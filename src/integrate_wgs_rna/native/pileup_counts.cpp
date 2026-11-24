// Python 생산기의 CIGAR 의미를 유지함. HTSlib pileup 검증기와 집계 루프를 공유하지 않음.
#include <htslib/sam.h>
#include <algorithm>
#include <climits>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>

struct Output {
    uint32_t *bases;
    uint64_t *mapq, *baseq;
    uint32_t *baseq_n, *forward, *deletion, *insertion;
    uint32_t *deletion_event, *max_deletion, *supplementary;
};

struct Context {
    samFile *file = nullptr;
    sam_hdr_t *header = nullptr;
    hts_idx_t *index = nullptr;
    hts_itr_t *iterator = nullptr;
    int tid = -1, min_mapq = 0, excluded_flags = 0;
    bool proper_pair = false;
    int64_t start = 0, end = 0;
    Output *output = nullptr;
    std::string error;
    ~Context() {
        if (iterator) hts_itr_destroy(iterator);
        if (index) hts_idx_destroy(index);
        if (header) sam_hdr_destroy(header);
        if (file) sam_close(file);
    }
};

static void message(char *buffer, size_t size, const std::string &value) {
    if (size) std::snprintf(buffer, size, "%s", value.c_str());
}

static bool increment(uint32_t &value, Context *context) {
    if (value == UINT32_MAX) {
        context->error = "uint32 observation overflow";
        return false;
    }
    ++value;
    return true;
}

extern "C" void *mfm_pileup_open(const char *path, const char *reference, const char *contig,
                                  int min_mapq, int excluded_flags, int proper_pair,
                                  char *error, size_t error_size) {
    try {
        auto context = std::make_unique<Context>();
        context->file = sam_open(path, "r");
        if (!context->file) throw std::runtime_error("alignment open failed");
        if (reference && *reference && hts_set_fai_filename(context->file, reference))
            throw std::runtime_error("reference setup failed");
        const int eof = hts_check_EOF(context->file);
        if (eof == 0 || eof == -1) throw std::runtime_error("alignment EOF check failed");
        context->header = sam_hdr_read(context->file);
        if (!context->header) throw std::runtime_error("alignment header failed");
        context->tid = sam_hdr_name2tid(context->header, contig);
        if (context->tid < 0) throw std::runtime_error("contig missing");
        context->index = sam_index_load(context->file, path);
        if (!context->index) throw std::runtime_error("alignment index missing");
        context->min_mapq = min_mapq;
        context->excluded_flags = excluded_flags;
        context->proper_pair = proper_pair;
        return context.release();
    } catch (const std::exception &exception) {
        message(error, error_size, exception.what());
        return nullptr;
    }
}

extern "C" int mfm_pileup_count(void *pointer, int64_t start, int64_t end, Output *out,
                                 char *error, size_t error_size) {
    auto *context = static_cast<Context *>(pointer);
    std::unique_ptr<bam1_t, decltype(&bam_destroy1)> record(bam_init1(), bam_destroy1);
    try {
        if (!context || !record || !out || start < 0 || end <= start)
            throw std::runtime_error("invalid pileup range");
        context->error.clear();
        if (context->iterator) hts_itr_destroy(context->iterator);
        context->iterator = sam_itr_queryi(context->index, context->tid, start, end);
        if (!context->iterator) throw std::runtime_error("alignment iterator creation failed");
        auto add = [&](uint32_t &value) {
            if (!increment(value, context)) throw std::runtime_error(context->error);
        };
        int status;
        while ((status = sam_itr_next(context->file, context->iterator, record.get())) >= 0) {
            const auto *read = record.get();
            const auto flags = read->core.flag;
            if ((flags & context->excluded_flags) || read->core.qual < context->min_mapq || !read->core.l_qseq)
                continue;
            if (context->proper_pair && (flags & BAM_FPAIRED) && !(flags & BAM_FPROPER_PAIR)) continue;
            const bool supplementary = flags & BAM_FSUPPLEMENTARY;
            const auto *cigar = bam_get_cigar(read);
            const auto *sequence = bam_get_seq(read);
            const auto *qualities = bam_get_qual(read);
            int64_t ref = read->core.pos, query = 0;
            for (uint32_t i = 0; i < read->core.n_cigar; ++i) {
                const auto op = bam_cigar_op(cigar[i]);
                const int64_t length = bam_cigar_oplen(cigar[i]);
                if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
                    const auto limit = std::min(end, ref + length);
                    for (int64_t position = std::max(start, ref); position < limit; ++position) {
                        const auto offset = position - start;
                        if (supplementary) { add(out->supplementary[offset]); continue; }
                        const auto q = query + position - ref;
                        if (q < 0 || q >= read->core.l_qseq) throw std::runtime_error("invalid query offset");
                        const int base = bam_seqi(sequence, q);
                        const int code = base == 1 ? 0 : base == 2 ? 1 : base == 4 ? 2 : base == 8 ? 3 : base == 15 ? 4 : -1;
                        if (code < 0) throw std::runtime_error("unsupported observed base");
                        add(out->bases[code*(end-start)+offset]);
                        out->mapq[offset] += read->core.qual;
                        if (!(flags & BAM_FREVERSE)) add(out->forward[offset]);
                        if (qualities[0] != 255 && qualities[q] != 255) {
                            out->baseq[offset] += qualities[q];
                            add(out->baseq_n[offset]);
                        }
                    }
                    ref += length; query += length;
                } else if (op == BAM_CINS || op == BAM_CDEL) {
                    const auto anchor = ref - 1;
                    if (!supplementary && ref > read->core.pos && anchor >= start && anchor < end) {
                        const auto offset = anchor - start;
                        if (op == BAM_CINS) add(out->insertion[offset]);
                        else {
                            add(out->deletion_event[offset]);
                            out->max_deletion[offset] = std::max(out->max_deletion[offset], static_cast<uint32_t>(length));
                        }
                    }
                    if (op == BAM_CINS) query += length;
                    else {
                        if (!supplementary) {
                            for (int64_t position = std::max(start,ref); position < std::min(end,ref+length); ++position)
                                add(out->deletion[position-start]);
                        }
                        ref += length;
                    }
                } else if (op == BAM_CREF_SKIP) ref += length;
                else if (op == BAM_CSOFT_CLIP) query += length;
                else if (op != BAM_CHARD_CLIP && op != BAM_CPAD)
                    throw std::runtime_error("unsupported CIGAR operation");
            }
        }
        if (status < -1) throw std::runtime_error("alignment iterator read failure");
        return 0;
    } catch (const std::exception &exception) {
        message(error, error_size, exception.what());
        return 1;
    }
}

extern "C" void mfm_pileup_close(void *pointer) { delete static_cast<Context *>(pointer); }
