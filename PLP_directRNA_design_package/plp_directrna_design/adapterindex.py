import logging
from typing import Optional, Tuple, Dict, List
import time

from cutadapt.adapters import Matchable, SingleAdapter, RemoveBeforeMatch, RemoveAfterMatch, PrefixAdapter, SuffixAdapter
from cutadapt.align import edit_environment, hamming_sphere


logger = logging.getLogger()

class AdapterIndex:
    """
    Index of multiple adapters

    Represent multiple adapters of the same type at once and use an index data structure
    to speed up matching. This is faster than iterating over multiple adapters.

    There are quite a few restrictions:
    - the error rate allows at most 2 mismatches
    - wildcards in the adapter are not allowed
    - wildcards in the read are not allowed

    Use the is_acceptable() method to check individual adapters.
    """

    AdapterIndexDict = Dict[str, Tuple[SingleAdapter, int, int]]

    def __init__(self, adapters, prefix: bool, keep_ambiguous: bool = False):
        """All given adapters must be of the same type"""
        if not adapters:
            raise ValueError("Adapter list is empty")
        for adapter in adapters:
            self._accept(adapter, prefix)
        self._adapters = adapters
        self._lengths, self._index, self._ambiguous = self._make_index(adapters, keep_ambiguous)
        logger.debug(
            "String lengths in the index: %s", sorted(self._lengths, reverse=True)
        )

        if len(self._lengths) == 1:
            self._length = self._lengths[0]
            self.match_to = self._match_to_one_length
        else:
            self.match_to = self._match_to_multiple_lengths
        if prefix:
            self._make_affix = self._make_prefix
            self._make_match = self._make_prefix_match
        else:
            self._make_affix = self._make_suffix
            self._make_match = self._make_suffix_match

    def __repr__(self):
        return f"{self.__class__.__name__}(adapters={self._adapters!r})"

    @staticmethod
    def _make_suffix(s, n):
        return s[-n:]

    @staticmethod
    def _make_prefix(s, n):
        return s[:n]

    @staticmethod
    def _make_prefix_match(adapter, length, score, errors, sequence):
        return RemoveBeforeMatch(
            astart=0,
            astop=len(adapter.sequence),
            rstart=0,
            rstop=length,
            score=score,
            errors=errors,
            adapter=adapter,
            sequence=sequence,
        )

    @staticmethod
    def _make_suffix_match(adapter, length, score, errors, sequence):
        return RemoveAfterMatch(
            astart=0,
            astop=len(adapter.sequence),
            rstart=len(sequence) - length,
            rstop=len(sequence),
            score=score,
            errors=errors,
            adapter=adapter,
            sequence=sequence,
        )

    @classmethod
    def _accept(cls, adapter: SingleAdapter, prefix: bool):
        """Raise a ValueError if the adapter is not acceptable"""
        if prefix and not isinstance(adapter, PrefixAdapter):
            raise ValueError("Only 5' anchored adapters are allowed")
        elif not prefix and not isinstance(adapter, SuffixAdapter):
            raise ValueError("Only 3' anchored adapters are allowed")
        if adapter.read_wildcards:
            raise ValueError("Wildcards in the read not supported")
        if adapter.adapter_wildcards:
            raise ValueError("Wildcards in the adapter not supported")
        k = int(len(adapter) * adapter.max_error_rate)
        if k > 3:
            raise ValueError("Error rate too high")

    @classmethod
    def is_acceptable(cls, adapter: SingleAdapter, prefix: bool):
        """
        Return whether this adapter is acceptable for being used in an index

        Adapters are not acceptable if they allow wildcards, allow too many errors,
        or would lead to a very large index.
        """
        try:
            cls._accept(adapter, prefix)
        except ValueError:
            return False
        return True

    @staticmethod
    def _make_index(adapters, keep_ambiguous: bool = False) -> Tuple[List[int], "AdapterIndexDict", int]:
        start_time = time.time()
        max_k = max(
            (
                int(adapter.max_error_rate * len(adapter.sequence))
                for adapter in adapters
                if adapter.indels
            ),
            default=0,
        )
        logger.info("Building index of %s adapters ...", len(adapters))
        if max_k == 3:
            logger.info(
                "Three errors and indels allowed for at least one of the adapter sequences: "
                "Indexing could take long and use a lot of memory. "
                "If this becomes a problem, try --no-indels and/or --no-index."
            )
        index: Dict[str, Tuple[SingleAdapter, int, int]] = dict()
        lengths = set()
        ambiguous = {}
        for adapter in adapters:
            sequence = adapter.sequence
            k = int(adapter.max_error_rate * len(sequence))

            if adapter.indels:
                for s, errors, matches in edit_environment(sequence, k):
                    if s in index:
                        other_adapter, other_errors, other_matches = index[s]
                        if matches < other_matches:
                            continue
                        if other_matches == matches and s not in ambiguous:
                            ambiguous[s] = (adapter, other_adapter, k, matches)
                    index[s] = (adapter, errors, matches)
                    lengths.add(len(s))
            else:
                n = len(sequence)
                for errors in range(k + 1):
                    matches = n - errors
                    for s in hamming_sphere(sequence, errors):
                        if s in index:
                            other_adapter, other_errors, other_matches = index[s]
                            if matches < other_matches:
                                continue
                            if other_matches == matches and s not in ambiguous:
                                ambiguous[s] = (adapter, other_adapter, k, matches)
                        index[s] = (adapter, errors, matches)
                lengths.add(n)

        if ambiguous and not keep_ambiguous:
            logger.warning(
                "WARNING: The adapters are too similar. When creating the index, "
                "%d ambiguous sequences were found that cannot be assigned uniquely.",
                len(ambiguous),
            )
            s = next(iter(ambiguous))
            adapter, other_adapter, k, matches = ambiguous[s]
            logger.warning(
                "WARNING: For example, %r, when found in a read, would result in "
                "%s matches for both %s %r and %s %r",
                s,
                matches,
                other_adapter.name,
                other_adapter.sequence,
                adapter.name,
                adapter.sequence,
            )
            logger.warning(
                "WARNING: Reads with ambiguous sequence will *not* be trimmed."
            )
            for s in ambiguous:
                del index[s]

        elapsed = time.time() - start_time
        logger.info("Built an index containing %s strings.", len(index))
        logger.debug("Building the index took %.1f s", elapsed)

        return sorted(lengths, reverse=True), index, len(ambiguous)

    def _match_to_one_length(self, sequence: str):
        """
        Match a query string against all adapters and return a Match that represents
        the best match or None if no match was found
        """
        affix = self._make_affix(sequence.upper(), self._length)
        if "N" in affix:
            result = self._lookup_with_n(affix)
            if result is None:
                return None
            adapter, e, m = result
        else:
            try:
                adapter, e, m = self._index[affix]
            except KeyError:
                return None
        return self._make_match(adapter, self._length, m, e, sequence)

    def _match_to_multiple_lengths(self, sequence: str):
        """
        Match the adapters against a string and return a Match that represents
        the best match or None if no match was found
        """
        affix = sequence.upper()

        # Check all the prefixes or suffixes (affixes) that could match
        best_adapter: Optional[SingleAdapter] = None
        best_length = 0
        best_m = -1
        best_e = 1000

        # Check successively shorter affixes
        for length in self._lengths:
            if length < best_m:
                # No chance of getting the same or a higher number of matches, so we can stop early
                break
            affix = self._make_affix(affix, length)
            if "N" in affix:
                result = self._lookup_with_n(affix)
                if result is None:
                    continue
                adapter, e, m = result
            else:
                try:
                    adapter, e, m = self._index[affix]
                except KeyError:
                    continue

            if m > best_m or (m == best_m and e < best_e):
                # TODO this could be made to work:
                # assert best_m == -1
                best_adapter = adapter
                best_e = e
                best_m = m
                best_length = length

        if best_m == -1:
            return None
        else:
            return self._make_match(best_adapter, best_length, best_m, best_e, sequence)

    def _lookup_with_n(self, affix):
        # N wildcards need to be counted as mismatches (read wildcards aren’t allowed).
        # We can thus look up an affix where we replace N with an arbitrary nucleotide.
        affix_without_n = affix.replace("N", "A")
        try:
            result = self._index[affix_without_n]
        except KeyError:
            return None

        # The looked up number of matches and errors is too low if
        # the adapter actually has an A where the N is in the query.
        # Fix this by re-doing the alignment.
        adapter = result[0]
        match = adapter.match_to(affix)
        if match is None:
            return None
        return adapter, match.errors, match.score