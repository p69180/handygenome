import itertools

import importlib

top_package_name = __name__.split(".")[0]
common = importlib.import_module(".".join([top_package_name, "common"]))


DEFAULT_FETCH_EXTEND_LENGTH = 10


def leftmost(vcfspec, fasta, fetch_extend_length=DEFAULT_FETCH_EXTEND_LENGTH):
    return indel_equivalents(vcfspec, fasta, fetch_extend_length)[0]


def rightmost(vcfspec, fasta, fetch_extend_length=DEFAULT_FETCH_EXTEND_LENGTH):
    return indel_equivalents(vcfspec, fasta, fetch_extend_length)[-1]


def check_equivalent(vcfspec1, vcfspec2, fasta):
    return leftmost(vcfspec1, fasta) == leftmost(vcfspec2, fasta)


def indel_equivalents(vcfspec, fasta, parsimoniously=True, fetch_extend_length=DEFAULT_FETCH_EXTEND_LENGTH):
    if parsimoniously:
        vcfspec = make_parsimonious(vcfspec)
    return shifting_equivalents(vcfspec=vcfspec, fasta=fasta, fetch_extend_length=fetch_extend_length)


def make_parsimonious(vcfspec):
    common.check_vcfspec_monoalt(vcfspec)

    pos = vcfspec.pos
    ref = vcfspec.ref
    alt = vcfspec.alts[0]

    # left
    min_length = min(len(ref), len(alt))
    action_idx = None
    for idx in range(0, min_length - 1):
        if ref[idx] == alt[idx]:
            action_idx = idx
            continue
        else:
            break

    if action_idx is not None:
        ref = ref[(action_idx + 1):]
        alt = alt[(action_idx + 1):]
        pos += (action_idx + 1)

    # right
    min_length = min(len(ref), len(alt))
    action_idx = None
    for idx in range(-1, -min_length, -1):
        if ref[idx] == alt[idx]:
            action_idx = idx
            continue
        else:
            break

    if action_idx is not None:
        ref = ref[:action_idx]
        alt = alt[:action_idx]

    # result
    return common.Vcfspec(vcfspec.chrom, pos, ref, (alt,))


def shifting_equivalents(vcfspec, fasta, fetch_extend_length):
    """Args:
        vcfspec: Must have single ALT.

    Returns:
        [vcfspec, vcfspec, ...]
            - each vcfspec is equivalent to the input vcfspec
            - each vcfspec is the most parsimonious form
            - sorted by pos (increasing order)

        Return value includes only input 'vcfspec' itself if one or more
            of the below is true:
            1) ref or alt sequence includes 'n' or 'N'
            2) input vcfspec is not a simple insertion or deletion (e.g. SNV, MNV, complex indel(delins))
    """
    ##### progressive fetch #####

    def progressive_fetch_forward(
        chrom,
        pos,
        fasta,
        fetch_extend_length,
        diff_seq,
        is_ins,
    ):
        extend_count = 0
        out_of_range = False
        chrom_length = fasta.lengths[fasta.references.index(chrom)]
        diff_seq_length = len(diff_seq)
        fetch_range_setter = (
            fetch_range_setter_ins_forward if is_ins else fetch_range_setter_del_forward
        )
        while True:
            if out_of_range:
                break

            extend_count += 1
            fetch_start, fetch_end = fetch_range_setter(
                pos, extend_count, fetch_extend_length, diff_seq_length
            )
            if fetch_end > chrom_length:
                out_of_range = True
                if fetch_start >= chrom_length:
                    continue

            result = fasta.fetch(chrom, fetch_start, fetch_end)
            yield result

    def progressive_fetch_backward(
        chrom,
        pos,
        fasta,
        fetch_extend_length,
    ):
        extend_count = 0
        out_of_range = False
        while True:
            if out_of_range:
                break

            extend_count += 1
            fetch_start, fetch_end = fetch_range_setter_backward(
                pos, extend_count, fetch_extend_length
            )
            if fetch_start < 0:
                out_of_range = True
                fetch_start = 0
                if fetch_end <= 0:
                    continue

            result = reversed(list(fasta.fetch(chrom, fetch_start, fetch_end)))
            yield result

    def fetch_range_setter_ins_forward(
        pos,
        extend_count,
        fetch_extend_length,
        diff_seq_length,
    ):
        fetch_start = pos + ((extend_count - 1) * fetch_extend_length)
        fetch_end = pos + (extend_count * fetch_extend_length)
        return fetch_start, fetch_end

    def fetch_range_setter_del_forward(
        pos,
        extend_count,
        fetch_extend_length,
        diff_seq_length,
    ):
        fetch_start = pos + diff_seq_length + ((extend_count - 1) * fetch_extend_length)
        fetch_end = pos + diff_seq_length + (extend_count * fetch_extend_length)
        return fetch_start, fetch_end

    def fetch_range_setter_backward(
        pos,
        extend_count,
        fetch_extend_length,
    ):
        fetch_start = pos - (extend_count * fetch_extend_length)
        fetch_end = pos - ((extend_count - 1) * fetch_extend_length)
        return fetch_start, fetch_end

    ##### navigate #####

    def navigate(
        chrom, pos, fasta, diff_seq, new_pos_list, fetch_extend_length, is_ins, forward
    ):
        """Args:
        is_ins: True for insertion, False for deletion
        forward: True for forward naviation (to 3' side on the
            positive strand), False for backward navigation
        """
        fetch_reservoir = list()
        fetched_seq = list()
        diff_seq_compared = diff_seq if forward else diff_seq[::-1]

        if forward:
            fetcher = progressive_fetch_forward(
                chrom,
                pos,
                fasta,
                fetch_extend_length,
                diff_seq,
                is_ins,
            )
            addnewpos = (
                navigate_helper_addnewpos_ins_forward
                if is_ins
                else navigate_helper_addnewpos_del_forward
            )
        else:
            fetcher = progressive_fetch_backward(
                chrom,
                pos,
                fasta,
                fetch_extend_length,
            )
            addnewpos = (
                navigate_helper_addnewpos_ins_backward
                if is_ins
                else navigate_helper_addnewpos_del_backward
            )

        for idx in itertools.count(0):
            navigate_helper_extend_fetch(fetch_reservoir, idx, diff_seq, fetcher)

            try:
                seq = fetch_reservoir[idx]
            except IndexError:
                break
            else:
                fetched_seq.append(seq)

            idx_remainder = idx % len(diff_seq)

            if fetched_seq[-1] != diff_seq_compared[idx_remainder]:
                break
            else:
                addnewpos(
                    idx,
                    idx_remainder,
                    fetch_reservoir,
                    pos,
                    diff_seq,
                    diff_seq_compared,
                    new_pos_list,
                    fetched_seq,
                )

    def navigate_helper_extend_fetch(fetch_reservoir, idx, diff_seq, fetcher):
        # "idx + 1 + len(diff_seq)" is the largest index of "fetch_reservoir", used in "navigate_helper_addnewpos_del_backward"
        limit_length = idx + 1 + len(diff_seq)
        while len(fetch_reservoir) <= limit_length:
            try:
                fetcher_output = next(fetcher)
            except StopIteration:
                break
            else:
                fetch_reservoir.extend(fetcher_output)

    def navigate_helper_addnewpos_ins_forward(
        idx,
        idx_remainder,
        fetch_reservoir,
        pos,
        diff_seq,
        diff_seq_compared,
        new_pos_list,
        fetched_seq,
    ):
        try:
            new_ref = fetch_reservoir[idx]
        except IndexError:
            pass
        else:
            new_pos = pos + (idx + 1)
            inserted_seq = (
                diff_seq_compared[idx_remainder + 1 :]
                + diff_seq_compared[: idx_remainder + 1]
            )
            new_alt = new_ref + "".join(inserted_seq)

            new_pos_list.append((new_pos, new_ref, new_alt))

    def navigate_helper_addnewpos_ins_backward(
        idx,
        idx_remainder,
        fetch_reservoir,
        pos,
        diff_seq,
        diff_seq_compared,
        new_pos_list,
        fetched_seq,
    ):
        try:
            new_ref = fetch_reservoir[idx + 1]
        except IndexError:
            pass
        else:
            new_pos = pos - (idx + 1)
            inserted_seq = list(
                reversed(
                    diff_seq_compared[idx_remainder + 1 :]
                    + diff_seq_compared[: idx_remainder + 1]
                )
            )
            new_alt = new_ref + "".join(inserted_seq)

            new_pos_list.append((new_pos, new_ref, new_alt))

    def navigate_helper_addnewpos_del_forward(
        idx,
        idx_remainder,
        fetch_reservoir,
        pos,
        diff_seq,
        diff_seq_compared,
        new_pos_list,
        fetched_seq,
    ):
        # overlapping deletion
        if idx < len(diff_seq):
            new_pos = pos + (idx + 1)
            deleted_seq = diff_seq_compared[idx + 1 :] + fetched_seq
            new_ref = diff_seq[idx] + "".join(deleted_seq)
            new_alt = new_ref[0]

            new_pos_list.append((new_pos, new_ref, new_alt))

        # non-overlapping deletion
        try:
            beyond_seq = get_beyond_seq(fetch_reservoir, idx, diff_seq)
        except IndexError:
            pass
        else:
            if nonoverlapping_deletion_check(
                diff_seq_compared, idx_remainder, beyond_seq
            ):
                new_pos = pos + len(diff_seq) + (idx + 1)
                deleted_seq = beyond_seq
                new_ref = fetch_reservoir[idx] + "".join(deleted_seq)
                # this line is not enclosed in try clause because the index used here is smaller than those used in "get_beyond_seq" function.
                new_alt = new_ref[0]

                new_pos_list.append((new_pos, new_ref, new_alt))

    def navigate_helper_addnewpos_del_backward(
        idx,
        idx_remainder,
        fetch_reservoir,
        pos,
        diff_seq,
        diff_seq_compared,
        new_pos_list,
        fetched_seq,
    ):
        # overlapping deletion
        if idx < len(diff_seq):
            try:
                seq = fetch_reservoir[idx + 1]
            except IndexError:
                pass
            else:
                new_pos = pos - (idx + 1)
                deleted_seq = list(reversed(diff_seq_compared[idx + 1 :] + fetched_seq))
                new_ref = seq + "".join(deleted_seq)
                new_alt = new_ref[0]

                new_pos_list.append((new_pos, new_ref, new_alt))

        # non-overlapping deletion
        try:
            beyond_seq = get_beyond_seq(fetch_reservoir, idx, diff_seq)
        except IndexError:
            pass
        else:
            if nonoverlapping_deletion_check(
                diff_seq_compared, idx_remainder, beyond_seq
            ):
                try:
                    seq = fetch_reservoir[idx + len(diff_seq) + 1]
                except IndexError:
                    pass
                else:
                    new_pos = pos - (idx + 1) - len(diff_seq)
                    deleted_seq = list(reversed(beyond_seq))
                    new_ref = seq + "".join(deleted_seq)
                    new_alt = new_ref[0]

                    new_pos_list.append((new_pos, new_ref, new_alt))

    def get_beyond_seq(fetch_reservoir, idx, diff_seq):
        if len(fetch_reservoir) < idx + 1 + len(diff_seq):
            raise IndexError
        else:
            return fetch_reservoir[(idx + 1) : (idx + 1 + len(diff_seq))]

    def nonoverlapping_deletion_check(diff_seq_compared, idx_remainder, beyond_seq):
        return (
            diff_seq_compared[idx_remainder + 1 :] == beyond_seq[: -(idx_remainder + 1)]
        ) and (
            diff_seq_compared[: idx_remainder + 1] == beyond_seq[-(idx_remainder + 1) :]
        )

    ##################################################

    def get_result(chrom, pos, ref, alt, fasta, diff_seq, fetch_extend_length, is_ins):
        # get new_pos_list
        new_pos_list = list()
        navigate(
            chrom,
            pos,
            fasta,
            diff_seq,
            new_pos_list,
            fetch_extend_length,
            is_ins,
            forward=True,
        )
        navigate(
            chrom,
            pos,
            fasta,
            diff_seq,
            new_pos_list,
            fetch_extend_length,
            is_ins,
            forward=False,
        )
        # make result
        result = list()
        result.append(common.Vcfspec(chrom, pos, ref, [alt]))
        for new_pos, new_ref, new_alt in new_pos_list:
            result.append(common.Vcfspec(chrom, new_pos, new_ref, [new_alt]))
        result.sort(key=(lambda x: x.pos))

        return result

    # main
    common.check_vcfspec_monoallele(vcfspec)

    chrom, pos, ref, alts = vcfspec.get_tuple()
    alt = alts[0]

    if "n" not in ref and "N" not in ref and "n" not in alt and "N" not in alt:
        # if set('nN').isdisjoint(ref) and set('nN').isdisjoint(alt):
        mttype = vcfspec.get_mttype_firstalt()
        if mttype in ("ins", "del"):
            is_ins = mttype == "ins"
            diff_seq = list(alt[1:]) if is_ins else list(ref[1:])
            result = get_result(
                chrom, pos, ref, alt, fasta, diff_seq, fetch_extend_length, is_ins
            )
        else:
            result = [vcfspec]
    else:
        result = [vcfspec]

    return result


"""
Below functions are kept for demonstration of the basic concept.

def _check_commutative_for_proof(refseq, insseq):
    def cond1(refseq, insseq):
        n = 0
        while True:
            n += 1
            end = n*len(insseq)
            start = (n-1)*len(insseq)
            subseq = refseq[start:end]
            if subseq != insseq:
                return False

            if end in range(len(refseq) - len(insseq) + 1, len(refseq) + 1):
                break

        return end

    def cond2(refseq, insseq):
        return refseq[ -len(insseq) : ] == insseq

    def cond3(refseq, insseq, end):
        #assert end < len(refseq)
        idx = len(refseq) - end
        return refseq[ -idx : ] == insseq[ : idx ]

    def main(refseq, insseq):
        end = cond1(refseq, insseq)
        if end == len(refseq):
            return True
        else:
            if cond2(refseq, insseq):
                return cond3(refseq, insseq, end)
            else:
                return False

    return main(refseq, insseq)
"""
