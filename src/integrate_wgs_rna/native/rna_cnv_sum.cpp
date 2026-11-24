// bedcov 기본 필터·D/N 포함 규칙의 구간합을 희소 endpoint에서 계산함.
#include <htslib/sam.h>
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
struct Region{int64_t begin,end;std::string gene;};
int main(int argc,char**argv) {
    samFile* in=nullptr;sam_hdr_t* header=nullptr;hts_idx_t* index=nullptr;hts_itr_t* iter=nullptr;bam1_t* read=nullptr;int rc=1;
    try {
        if(argc!=6)throw std::runtime_error("bam reference|- contig bed decode_threads");
        int threads=std::stoi(argv[5]);if(threads<0)throw std::runtime_error("invalid decoder count");
        std::ifstream bed(argv[4]);if(!bed)throw std::runtime_error("BED missing");
        std::vector<Region> rows;std::vector<int64_t> points;std::string line;
        while(std::getline(bed,line)) {
            Region row;std::vector<std::string> fields;size_t begin=0;
            for(size_t end;(end=line.find('\t',begin))!=std::string::npos;begin=end+1)fields.push_back(line.substr(begin,end-begin));
            fields.push_back(line.substr(begin));
            if(fields.size()!=4||fields[0]!=argv[3])throw std::runtime_error("invalid BED row");
            row.begin=std::stoll(fields[1]);row.end=std::stoll(fields[2]);row.gene=fields[3];
            if(row.begin<0||row.end<row.begin)throw std::runtime_error("invalid BED coordinate");
            points.push_back(row.begin);points.push_back(row.end);rows.push_back(std::move(row));
        }
        if(!bed.eof()||rows.empty())throw std::runtime_error("empty/corrupt BED");
        std::sort(points.begin(),points.end());points.erase(std::unique(points.begin(),points.end()),points.end());
        std::vector<int64_t> counts(points.size()+1,0),values(points.size(),0);
        std::vector<__int128> offsets(points.size()+1,0);
        in=sam_open(argv[1],"r");if(!in)throw std::runtime_error("alignment open failed");
        if(std::string(argv[2])!="-"&&hts_set_fai_filename(in,argv[2])!=0)throw std::runtime_error("reference setup failed");
        int eof=hts_check_EOF(in);if(eof==0||eof==-1)throw std::runtime_error("truncated alignment EOF");
        if(threads&&hts_set_threads(in,threads)!=0)throw std::runtime_error("decoder setup failed");
        const char* selected=std::getenv("RNA_CNV_CRAM_FIELDS");std::string fields=selected?selected:"coverage";
        if(fields!="full"&&fields!="coverage")throw std::runtime_error("RNA_CNV_CRAM_FIELDS must be full/coverage");
        if(fields=="coverage"&&hts_get_format(in)->format==cram) {
            // SEQ는 reference MD5 검증을 보존하기 위해 남기고 QUAL/aux/mate/name만 생략함.
            int required=SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_SEQ;
            if(hts_set_opt(in,CRAM_OPT_REQUIRED_FIELDS,required)!=0)throw std::runtime_error("CRAM required fields failed");
        }
        header=sam_hdr_read(in);if(!header)throw std::runtime_error("alignment header failed");
        int tid=sam_hdr_name2tid(header,argv[3]);if(tid<0)throw std::runtime_error("BED contig missing");
        index=sam_index_load(in,argv[1]);if(!index)throw std::runtime_error("alignment index missing");
        iter=sam_itr_queryi(index,tid,points.front(),points.back());if(!iter)throw std::runtime_error("iterator failed");
        read=bam_init1();if(!read)throw std::bad_alloc();int status;uint64_t fetched=0,retained=0;
        while((status=sam_itr_next(in,iter,read))>=0) {
            ++fetched;if(read->core.flag&1796)continue;
            int64_t span=bam_cigar2rlen(read->core.n_cigar,bam_get_cigar(read));
            if(span<=0)continue;
            for(uint32_t i=0;i<read->core.n_cigar;++i)if(bam_cigar_op(bam_get_cigar(read)[i])==BAM_CBACK)throw std::runtime_error("unsupported backward CIGAR");
            // shortcut: contig당 retained INT32_MAX 이하는 bedcov depth cap 미작동을 보장함. 초심도 입력은 cap 정책 검증 후 확장함.
            if(++retained>uint64_t(INT32_MAX))throw std::runtime_error("possible bedcov depth-cap overflow");
            int64_t start=read->core.pos,end=start+span;
            size_t first=std::lower_bound(points.begin(),points.end(),start)-points.begin();
            size_t last=std::lower_bound(points.begin(),points.end(),end)-points.begin();
            ++counts[first];offsets[first]-=start;--counts[last];offsets[last]+=end;
        }
        if(status<-1)throw std::runtime_error("truncated/corrupt alignment body");
        int64_t slope=0;__int128 intercept=0;
        for(size_t i=0;i<points.size();++i) {
            slope+=counts[i];intercept+=offsets[i];__int128 value=__int128(slope)*points[i]+intercept;
            if(value<0||value>INT64_MAX)throw std::runtime_error("coverage sum overflow");values[i]=int64_t(value);
        }
        for(const Region& row:rows) {
            auto lo=std::lower_bound(points.begin(),points.end(),row.begin)-points.begin();
            auto hi=std::lower_bound(points.begin(),points.end(),row.end)-points.begin();
            int64_t sum=values[hi]-values[lo];if(sum<0)throw std::runtime_error("negative coverage");
            std::cout<<argv[3]<<'\t'<<row.begin<<'\t'<<row.end<<'\t'<<row.gene<<'\t'<<sum<<'\n';
        }
        if(!std::cout)throw std::runtime_error("coverage output failed");
        std::cerr<<"CNV_SUM records="<<fetched<<" retained="<<retained<<" endpoints="<<points.size()<<" fields="<<fields<<'\n';rc=0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';}
    if(read)bam_destroy1(read);if(iter)hts_itr_destroy(iter);if(index)hts_idx_destroy(index);if(header)sam_hdr_destroy(header);if(in&&sam_close(in)<0)rc=1;return rc;
}
