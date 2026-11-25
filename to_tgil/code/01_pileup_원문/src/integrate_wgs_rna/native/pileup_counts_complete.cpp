// BAM/CRAM 한 번 해독하고 CIGAR 한 번 순회하며 공통·HP1·HP2를 집계함.
#include "pileup_complete.hpp"
using namespace complete_pileup;

extern "C" void *mfm_complete_open(const char *path,const char *reference,const char *contig,
                                    int mq,int flags,int proper,char *error,size_t size) {
    try{return open(path,reference,contig,mq,flags,proper);}
    catch(const std::exception &exception){message(error,size,exception.what());return nullptr;}
}
extern "C" int mfm_complete_count(void *pointer,int64_t start,int64_t end,Output *output,char *error,size_t size) {
    auto *context=static_cast<Context *>(pointer);
    try {
        if(!context)throw std::runtime_error("closed count context");
        context->range(start,end,output);
        std::unique_ptr<bam1_t,decltype(&bam_destroy1)> read(bam_init1(),bam_destroy1);
        if(!read)throw std::bad_alloc();
        int status;
        while((status=sam_itr_next(context->file,context->iterator,read.get()))>=0)
            if(context->accepted(read.get()))context->cigar(read.get(),false);
        if(status< -1)throw std::runtime_error("alignment iterator read failure");
        return 0;
    } catch(const std::exception &exception){message(error,size,exception.what());return 1;}
}
extern "C" void mfm_complete_close(void *pointer){delete static_cast<Context *>(pointer);}
