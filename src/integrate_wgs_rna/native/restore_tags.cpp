// 정렬 순서와 모든 read 필드를 보존하며 기존 RG/OD 복원만 수행함.
#include <htslib/sam.h>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <cstring>
int main(int argc,char**argv) {
  if(argc!=8){std::cerr<<"input output header rg mode output_format decode_threads\n";return 2;}
  samFile *in=nullptr,*out=nullptr;sam_hdr_t *old=nullptr,*header=nullptr;bam1_t *b=nullptr;int result=1;
  try {
    in=sam_open(argv[1],"rb");if(!in)throw std::runtime_error("input open failed");
    const int threads=std::stoi(argv[7]);if(threads<0)throw std::runtime_error("invalid decode threads");
    if(threads&&hts_set_threads(in,threads)<0)throw std::runtime_error("decode pool failed");
    old=sam_hdr_read(in);if(!old)throw std::runtime_error("input header failed");
    std::ifstream f(argv[3]);std::string text((std::istreambuf_iterator<char>(f)),{});
    header=sam_hdr_parse(text.size(),text.c_str());if(!header)throw std::runtime_error("target header failed");
    out=sam_open(argv[2],argv[6]);if(!out||sam_hdr_write(out,header)<0)throw std::runtime_error("output header failed");
    b=bam_init1();if(!b)throw std::bad_alloc();int status;uint64_t records=0;
    while((status=sam_read1(in,old,b))>=0) {
      uint8_t *rg=bam_aux_get(b,"RG");if(rg&&bam_aux_del(b,rg)<0)throw std::runtime_error("RG remove failed");
      if(bam_aux_append(b,"RG",'Z',std::strlen(argv[4])+1,(const uint8_t*)argv[4])<0)throw std::runtime_error("RG append failed");
      if(std::strcmp(argv[5],"unset")) {
        uint8_t *od=bam_aux_get(b,"OD");
        if(od){
          if(!std::strchr("cCsSiI",*od))throw std::runtime_error("OD must be integer");
          if(bam_aux2i(od)==1&&!std::strcmp(argv[5],"preserve_source"))b->core.flag|=BAM_FDUP;
          if(bam_aux_del(b,od)<0)throw std::runtime_error("OD remove failed");
        }
      }
      if(sam_write1(out,header,b)<0)throw std::runtime_error("record write failed");
      ++records;
    }
    if(status<-1)throw std::runtime_error("truncated or corrupt input");
    if(sam_close(out)<0){out=nullptr;throw std::runtime_error("output close failed");}out=nullptr;result=0;
    std::cerr<<"MFM_RESTORE_RECORDS="<<records<<'\n';
  }catch(const std::exception&e){std::cerr<<e.what()<<'\n';}
  if(b)bam_destroy1(b);if(old)sam_hdr_destroy(old);if(header)sam_hdr_destroy(header);
  if(in&&sam_close(in)<0)result=1;if(out)sam_close(out);return result;
}
