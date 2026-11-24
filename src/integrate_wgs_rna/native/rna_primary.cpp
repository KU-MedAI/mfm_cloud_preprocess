// 독립 전수 flag QC. primary mapped만 세며 QC-fail은 포함함.
#include <htslib/sam.h>
#include <cstdint>
#include <iostream>
int main(int argc,char**argv) {
    if(argc!=3)return 2;samFile* in=sam_open(argv[1],"r");if(!in)return 1;
    int threads=std::stoi(argv[2]);if(threads && hts_set_threads(in,threads)<0){sam_close(in);return 1;}
    auto* h=sam_hdr_read(in);auto* b=bam_init1();if(!h || !b){if(h)sam_hdr_destroy(h);if(b)bam_destroy1(b);sam_close(in);return 1;}
    uint64_t primary=0,duplicates=0;int status;
    while((status=sam_read1(in,h,b))>=0)if(!(b->core.flag&2308)){++primary;duplicates+=bool(b->core.flag&BAM_FDUP);}
    bam_destroy1(b);sam_hdr_destroy(h);int closed=sam_close(in);if(status<-1 || closed<0)return 1;
    std::cout<<"["<<duplicates<<","<<primary<<"]\n";return 0;
}
