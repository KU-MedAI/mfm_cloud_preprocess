// Primary base/deletion span은 생산 CIGAR 루프와 별개인 HTSlib pileup으로 검증함.
#include "pileup_complete.hpp"
using namespace complete_pileup;

static int next_read(void *pointer,bam1_t *read) {
    auto *context=static_cast<Context *>(pointer);
    try {
        int status;
        while((status=sam_itr_next(context->file,context->iterator,read))>=0) {
            if(!context->accepted(read))continue;
            // 단일 signed indel로는 같은 anchor의 I+D를 보존할 수 없음.
            context->cigar(read,true);
            return status;
        }
        if(status< -1)context->error="alignment iterator read failure";
        return status;
    } catch(const std::exception &exception){context->error=exception.what();return -2;}
}
static int initialize_phase(void *,const bam1_t *read,bam_pileup_cd *data) {
    data->i=Context::hp(read);return 0;
}
extern "C" void *mfm_complete_open(const char *path,const char *reference,const char *contig,
                                    int mq,int flags,int proper,char *error,size_t size) {
    try{return open(path,reference,contig,mq,flags,proper);}
    catch(const std::exception &exception){message(error,size,exception.what());return nullptr;}
}
extern "C" int mfm_complete_count(void *pointer,int64_t start,int64_t end,Output *output,char *error,size_t size) {
    auto *context=static_cast<Context *>(pointer);bam_plp_t pileup=nullptr;
    try {
        if(!context)throw std::runtime_error("closed count context");
        context->range(start,end,output);pileup=bam_plp_init(next_read,context);
        if(!pileup)throw std::bad_alloc();
        bam_plp_set_maxcnt(pileup,INT_MAX);bam_plp_constructor(pileup,initialize_phase);
        int tid,n;hts_pos_t position;const bam_pileup1_t *column;
        while((column=bam_plp64_auto(pileup,&tid,&position,&n))) {
            if(position<start)continue;if(position>=end)break;
            if(n==INT_MAX)throw std::runtime_error("integer depth limit reached; refusing silent cap");
            for(int i=0;i<n;++i) {
                const auto &entry=column[i];if(entry.is_refskip)continue;
                int phase=int(entry.cd.i);
                if(entry.is_del)context->deletion(entry.b,position,phase);
                else context->base(entry.b,position,entry.qpos,phase);
            }
        }
        if(!context->error.empty())throw std::runtime_error(context->error);
        bam_plp_destroy(pileup);return 0;
    } catch(const std::exception &exception) {
        if(pileup)bam_plp_destroy(pileup);message(error,size,exception.what());return 1;
    }
}
extern "C" void mfm_complete_close(void *pointer){delete static_cast<Context *>(pointer);}
