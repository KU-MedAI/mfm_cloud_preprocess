// shard 해석과 기존 fragment 배정을 통합하고 유전자별 작은 count만 Python에 반환함.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <ctime>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace ase_fused {
using Clock=std::chrono::steady_clock;
static double seconds(Clock::time_point t){return std::chrono::duration<double>(Clock::now()-t).count();}
static double cpu(){return double(std::clock())/CLOCKS_PER_SEC;}
static int64_t integer(PyObject* p){int64_t x=PyLong_AsLongLong(p);if(PyErr_Occurred())throw std::runtime_error("정수 범위 오류");return x;}
static std::string text(PyObject* p){Py_ssize_t n;const char* s=PyUnicode_AsUTF8AndSize(p,&n);if(!s)throw std::runtime_error("문자열 오류");return std::string(s,size_t(n));}
static void utf8(const char* p,size_t n){
    for(size_t i=0;i<n;){unsigned a=(unsigned char)p[i++];if(a<128)continue;
        unsigned lo=128,hi=191,count;
        if(a>=194&&a<=223)count=1;
        else if(a>=224&&a<=239){count=2;if(a==224)lo=160;if(a==237)hi=159;}
        else if(a>=240&&a<=244){count=3;if(a==240)lo=144;if(a==244)hi=143;}
        else throw std::runtime_error("잘못된 UTF-8 ASE 이름");
        if(count>n-i)throw std::runtime_error("잘린 UTF-8 ASE 이름");
        unsigned b=(unsigned char)p[i++];if(b<lo||b>hi)throw std::runtime_error("잘못된 UTF-8 ASE 이름");
        while(--count){b=(unsigned char)p[i++];if(b<128||b>191)throw std::runtime_error("잘못된 UTF-8 ASE 이름");}
    }
}
struct Block{int64_t begin,end;};
struct Mate{bool present=false;int contig=0;__int128 score=0;std::vector<Block> blocks;};
using Slot=std::array<Mate,2>;
using Table=std::unordered_map<std::string,Slot>;
struct Exon{int64_t begin,end;int gene;};
using Bins=std::unordered_map<int,std::unordered_map<int64_t,std::vector<Exon>>>;
struct Evidence{std::vector<int64_t> starts,ends;int64_t maxlen=0;};
struct Engine{
    std::unordered_map<std::string,int> contig_ids,gene_ids;
    std::vector<std::string> genes;
    Bins exons[2];std::unordered_map<int,Evidence> evidence[2][2];
    int mapq=255;bool duplicates=false,both=false;double scale=0;PyObject* split=nullptr;
    ~Engine(){Py_XDECREF(split);}
    int contig(const std::string& s){auto it=contig_ids.find(s);if(it!=contig_ids.end())return it->second;int id=int(contig_ids.size());contig_ids.emplace(s,id);return id;}
    int gene(const std::string& s){auto it=gene_ids.find(s);if(it!=gene_ids.end())return it->second;int id=int(genes.size());genes.push_back(s);gene_ids.emplace(s,id);return id;}
};
static void capsule_free(PyObject* capsule){auto* e=static_cast<Engine*>(PyCapsule_GetPointer(capsule,"mfm_ase_engine"));if(e)delete e;else PyErr_Clear();}
static void load_exons(Engine& e,PyObject* source,int hap){
    if(!PyDict_Check(source))throw std::runtime_error("ASE exon index가 dict가 아님");
    Py_ssize_t pos=0;PyObject *key,*value;
    while(PyDict_Next(source,&pos,&key,&value)){
        if(!PyTuple_Check(key)||PyTuple_Size(key)!=2)throw std::runtime_error("ASE exon key 오류");
        int c=e.contig(text(PyTuple_GET_ITEM(key,0)));int64_t bin=integer(PyTuple_GET_ITEM(key,1));
        Object rows(checked(PySequence_Fast(value,"exon 행 오류")));
        auto& dest=e.exons[hap][c][bin];dest.reserve(PySequence_Fast_GET_SIZE(rows.value));
        for(Py_ssize_t i=0;i<PySequence_Fast_GET_SIZE(rows.value);++i){
            Object row(checked(PySequence_Fast(PySequence_Fast_GET_ITEM(rows.value,i),"exon 행 오류")));
            if(PySequence_Fast_GET_SIZE(row.value)!=3)throw std::runtime_error("exon 열 오류");
            dest.push_back({integer(PySequence_Fast_GET_ITEM(row.value,0)),integer(PySequence_Fast_GET_ITEM(row.value,1)),e.gene(text(PySequence_Fast_GET_ITEM(row.value,2)))});
        }
    }
}
static void load_evidence(Engine& e,PyObject* source){
    if(!PyDict_Check(source))throw std::runtime_error("ASE evidence index가 dict가 아님");
    Py_ssize_t pos=0;PyObject *key,*value;
    while(PyDict_Next(source,&pos,&key,&value)){
        if(!PyTuple_Check(key)||PyTuple_Size(key)!=3)throw std::runtime_error("ASE evidence key 오류");
        std::string kind=text(PyTuple_GET_ITEM(key,0));int hap=int(integer(PyTuple_GET_ITEM(key,1)))-1;
        if((kind!="good"&&kind!="bad")||hap<0||hap>1)throw std::runtime_error("ASE evidence 종류 오류");
        int c=e.contig(text(PyTuple_GET_ITEM(key,2)));Object row(checked(PySequence_Fast(value,"evidence 행 오류")));
        if(PySequence_Fast_GET_SIZE(row.value)!=3)throw std::runtime_error("ASE evidence 열 오류");
        Object starts(checked(PySequence_Fast(PySequence_Fast_GET_ITEM(row.value,0),"evidence start 오류")));
        Object ends(checked(PySequence_Fast(PySequence_Fast_GET_ITEM(row.value,1),"evidence end 오류")));
        if(PySequence_Fast_GET_SIZE(starts.value)!=PySequence_Fast_GET_SIZE(ends.value))throw std::runtime_error("ASE evidence 길이 오류");
        Evidence item;item.maxlen=integer(PySequence_Fast_GET_ITEM(row.value,2));
        for(Py_ssize_t i=0;i<PySequence_Fast_GET_SIZE(starts.value);++i){item.starts.push_back(integer(PySequence_Fast_GET_ITEM(starts.value,i)));item.ends.push_back(integer(PySequence_Fast_GET_ITEM(ends.value,i)));}
        e.evidence[kind=="bad"][hap][c]=std::move(item);
    }
}
static PyObject* prepare(PyObject*,PyObject* args){
    PyObject *ex1,*ex2,*evidence,*split;int both,mapq,duplicates;double scale;
    if(!PyArg_ParseTuple(args,"OOOpipdO",&ex1,&ex2,&evidence,&both,&mapq,&duplicates,&scale,&split))return nullptr;
    Engine* e=nullptr;
    try{
        if(!PyCallable_Check(split))throw std::runtime_error("split callback이 함수가 아님");
        if(!std::isfinite(scale)||std::abs(scale)>1e12){PyErr_SetString(PyExc_OverflowError,"native ASE 길이 scale 범위 초과");return nullptr;}
        e=new Engine;e->both=both;e->mapq=mapq;e->duplicates=duplicates;e->scale=scale;Py_INCREF(split);e->split=split;
        load_exons(*e,ex1,0);load_exons(*e,ex2,1);load_evidence(*e,evidence);
        PyObject* result=PyCapsule_New(e,"mfm_ase_engine",capsule_free);if(!result)delete e;return result;
    }catch(const std::exception& error){delete e;if(!PyErr_Occurred())PyErr_SetString(PyExc_ValueError,error.what());return nullptr;}
}
static void read_file(PyObject* rawpath,Engine& e,int hap,std::array<Table,4>& tables,std::unordered_set<std::string>& dup,uint64_t& dupcount){
    Object path;if(!PyUnicode_FSConverter(rawpath,&path.value))throw std::runtime_error("shard 경로 오류");
    if(std::strlen(PyBytes_AS_STRING(path.value))!=size_t(PyBytes_GET_SIZE(path.value)))throw std::runtime_error("경로의 NUL 오류");
    Input input;input.fd=open(PyBytes_AS_STRING(path.value),O_RDONLY);
    if(input.fd<0){PyErr_SetFromErrnoWithFilenameObject(PyExc_OSError,rawpath);throw std::runtime_error("shard 열기 실패");}
    struct stat st;if(fstat(input.fd,&st)<0){PyErr_SetFromErrno(PyExc_OSError);throw std::runtime_error("shard stat 실패");}
    if(st.st_size<12)throw std::runtime_error("잘린 ASE shard");input.size=size_t(st.st_size);
    void* mapping=mmap(nullptr,input.size,PROT_READ,MAP_PRIVATE,input.fd,0);
    if(mapping==MAP_FAILED){PyErr_SetFromErrno(PyExc_OSError);throw std::runtime_error("shard mmap 실패");}
    input.data=static_cast<const char*>(mapping);
    if(std::memcmp(input.take(8),"MFMASE1\0",8))throw std::runtime_error("ASE shard magic 오류");
    uint32_t n=uint32_t(input.number(4));if(n>1000000)throw std::runtime_error("ASE reference 개수 오류");std::vector<int> refs;refs.reserve(n);
    for(uint32_t i=0;i<n;++i){size_t len=size_t(input.number(4));if(len>1<<20)throw std::runtime_error("ASE reference 길이 오류");const char* raw=input.take(len);utf8(raw,len);refs.push_back(e.contig(std::string(raw,len)));}
    while(input.pos<input.size){
        uint32_t len=uint32_t(input.number(4));uint16_t flag=uint16_t(input.number(2));int mq=int(input.number(1)),has=int(input.number(1));
        int64_t score=signed64(input.number(8));int32_t tid=int32_t(input.number(4));uint32_t blocks=uint32_t(input.number(4));
        if(len>1<<20||blocks>1<<20||tid<0||uint32_t(tid)>=n)throw std::runtime_error("ASE record 범위 오류");
        const char* rawname=input.take(len);utf8(rawname,len);std::string name(rawname,len);
        const auto* raw=reinterpret_cast<const unsigned char*>(input.take(size_t(blocks)*16));
        if(flag&1024){++dupcount;dup.insert(name);if(!e.duplicates)continue;}
        if(!has||(mq==0&&mq<e.mapq))continue;
        auto& slot=tables[hap*2+int(mq<e.mapq)].try_emplace(std::move(name)).first->second;
        auto& m=slot[(flag&64)?0:1];if(m.present&&score<=m.score)continue;
        m.present=true;m.score=score;m.contig=refs[tid];m.blocks.clear();m.blocks.reserve(blocks);
        for(uint32_t i=0;i<blocks;++i){uint64_t a=0,b=0;for(int k=0;k<8;++k){a|=uint64_t(raw[i*16+k])<<(8*k);b|=uint64_t(raw[i*16+8+k])<<(8*k);}m.blocks.push_back({signed64(a),signed64(b)});}
    }
}
static uint64_t correct(Table& table,double scale){
    if(scale==0)return 0;uint64_t changed=0;
    for(auto& row:table){auto& a=row.second[0];auto& b=row.second[1];if(a.present&&b.present&&a.contig!=b.contig)continue;
        const auto& x=a.blocks;const auto& y=b.blocks;int64_t start,end;
        if(!x.empty()&&!y.empty()){start=std::min(x.front().begin,y.front().begin);end=std::max(x.back().end,y.back().end);}
        else if(!x.empty()){start=x.front().begin;end=x.back().end;}
        else if(!y.empty()){start=y.front().begin;end=y.back().end;}
        else throw std::runtime_error("min() arg is an empty sequence");
        __int128 span=__int128(end)-start;double term=span<=0?0:std::ceil(std::log2(double(span))*scale-0.5);
        if(!std::isfinite(term)||term<=-std::ldexp(1.0,62)||term>=std::ldexp(1.0,62)){PyErr_SetString(PyExc_OverflowError,"native ASE 길이 항 범위 초과");throw std::runtime_error("길이 항 범위 오류");}
        if(term){if(a.present)a.score-=int64_t(term);if(b.present)b.score-=int64_t(term);++changed;}
    }return changed;
}
static bool hit(const Engine& e,bool bad,int hap,const Slot& slot){
    for(const auto& m:slot)if(m.present){auto it=e.evidence[bad][hap].find(m.contig);if(it==e.evidence[bad][hap].end())continue;const auto& v=it->second;
        for(const auto& b:m.blocks){auto first=std::lower_bound(v.starts.begin(),v.starts.end(),__int128(b.begin)-v.maxlen,[](int64_t x,__int128 y){return x<y;});auto last=std::lower_bound(v.starts.begin(),v.starts.end(),b.end);
            for(auto p=first;p<last;++p)if(v.ends[size_t(p-v.starts.begin())]>b.begin)return true;
        }
    }return false;
}
static int64_t bin(int64_t x){return x/16384-(x%16384<0);}
static std::unordered_set<int> genes_of(const Engine& e,int hap,const Slot& slot){
    std::unordered_set<int> genes;
    for(const auto& m:slot)if(m.present){auto c=e.exons[hap].find(m.contig);if(c==e.exons[hap].end())continue;
        for(const auto& b:m.blocks){if(b.end<=b.begin)continue;int64_t last=bin(b.end-1);
            for(int64_t k=bin(b.begin);;){auto found=c->second.find(k);if(found!=c->second.end())for(const auto& exon:found->second)if(exon.begin<b.end&&exon.end>b.begin)genes.insert(exon.gene);if(k==last)break;++k;}
        }
    }return genes;
}
struct Counts{std::array<std::unordered_map<int,uint64_t>,8> values;std::array<uint64_t,8> stats{};};
static void assign_one(const Engine& e,const std::string& name,const Slot* a,const Slot* b,Counts& result){
    int win;bool paired=a&&b;const Slot* slot;
    if(paired){int common=0;__int128 s1=0,s2=0;for(int i=0;i<2;++i)if((*a)[i].present&&(*b)[i].present){++common;s1+=(*a)[i].score;s2+=(*b)[i].score;}
        if(!common){++result.stats[4];return;}if(s1==s2){++result.stats[0];return;}win=s1>s2?0:1;slot=win==0?a:b;
        if(hit(e,true,0,*a)||hit(e,true,1,*b)){++result.stats[6];return;}
    }else{if(e.both){++result.stats[5];return;}win=a?0:1;slot=a?a:b;if(hit(e,true,win,*slot)){++result.stats[6];return;}}
    if(!hit(e,false,win,*slot)){++result.stats[7];return;}
    Object key(checked(PyUnicode_DecodeUTF8(name.data(),Py_ssize_t(name.size()),"strict")));
    Object bucket(checked(PyObject_CallOneArg(e.split,key.value)));long half=PyLong_AsLong(bucket.value);
    if(PyErr_Occurred())throw std::runtime_error("split 오류");if(half<0||half>1)throw std::runtime_error("split 범위 오류");
    if(!paired)++result.stats[win+1];
    for(int gene:genes_of(e,win,*slot)){++result.values[win][gene];++result.values[4+win*2+half][gene];if(paired)++result.values[2+win][gene];}
    ++result.stats[3];
}
static Counts assign(const Engine& e,const Table& a,const Table& b){
    Counts out;for(const auto& row:a){auto it=b.find(row.first);assign_one(e,row.first,&row.second,it==b.end()?nullptr:&it->second,out);}
    for(const auto& row:b)if(a.find(row.first)==a.end())assign_one(e,row.first,nullptr,&row.second,out);return out;
}
static PyObject* dict_counts(const Engine& e,const Counts& c){
    Object result(checked(PyDict_New()));std::array<PyObject*,8> values{};std::vector<Object*> owned;
    try{for(int i=0;i<8;++i){auto* v=new Object(checked(PyDict_New()));owned.push_back(v);values[i]=v->value;
            for(const auto& row:c.values[i]){const auto& name=e.genes[row.first];Object key(checked(PyUnicode_DecodeUTF8(name.data(),Py_ssize_t(name.size()),"strict")));Object n(checked(PyLong_FromUnsignedLongLong(row.second)));success(PyDict_SetItem(v->value,key.value,n.value));}}
        const char* names[]={"a","b","both_a","both_b"};for(int i=0;i<4;++i)success(PyDict_SetItemString(result.value,names[i],values[i]));
        Object sa(checked(PyTuple_Pack(2,values[4],values[5]))),sb(checked(PyTuple_Pack(2,values[6],values[7])));success(PyDict_SetItemString(result.value,"split_a",sa.value));success(PyDict_SetItemString(result.value,"split_b",sb.value));
    }catch(...){for(auto* p:owned)delete p;throw;}
    for(auto* p:owned)delete p;return result.release();
}
static PyObject* dict_stats(const Counts& c){
    const char* names[]={"tie","only1","only2","assigned","nocommon","dropped_single","drop_unvalidated_variant","drop_no_validated_variant"};Object out(checked(PyDict_New()));
    for(int i=0;i<8;++i){Object n(checked(PyLong_FromUnsignedLongLong(c.stats[i])));success(PyDict_SetItemString(out.value,names[i],n.value));}return out.release();
}
static PyObject* run(PyObject*,PyObject* args){
    PyObject *path1,*path2,*capsule;if(!PyArg_ParseTuple(args,"OOO",&path1,&path2,&capsule))return nullptr;
    auto* e=static_cast<Engine*>(PyCapsule_GetPointer(capsule,"mfm_ase_engine"));if(!e)return nullptr;
    try{
        auto started=Clock::now();double c0=cpu();std::array<Table,4> tables;std::unordered_set<std::string> d1,d2;uint64_t n1=0,n2=0;
        read_file(path1,*e,0,tables,d1,n1);read_file(path2,*e,1,tables,d2,n2);double read_seconds=seconds(started),read_cpu=cpu()-c0;
        auto filtered=Clock::now();double cf=cpu();std::unordered_set<std::string> united=d1;united.insert(d2.begin(),d2.end());uint64_t removed=0;
        if(!e->duplicates)for(auto& table:tables)for(const auto& name:united)removed+=table.erase(name);
        uint64_t corrected=0;for(auto& table:tables)corrected+=correct(table,e->scale);
        uint64_t total=tables[0].size();for(const auto& row:tables[2])if(tables[0].find(row.first)==tables[0].end())++total;
        uint64_t intersection=d1.size()+d2.size()-united.size();double filter_seconds=seconds(filtered),filter_cpu=cpu()-cf;
        auto assigned=Clock::now();double ca=cpu();Counts one=assign(*e,tables[0],tables[2]),mm=assign(*e,tables[1],tables[3]);
        Object a(checked(dict_counts(*e,one))),b(checked(dict_stats(one))),c(checked(dict_counts(*e,mm))),d(checked(dict_stats(mm)));
        Object audit(checked(Py_BuildValue("{s:K,s:K,s:K,s:K,s:K,s:K,s:K}","dup1",(unsigned long long)n1,"dup2",(unsigned long long)n2,"dup_union",(unsigned long long)united.size(),"dup_inconsistent",(unsigned long long)(d1.size()+d2.size()-2*intersection),"removed",(unsigned long long)removed,"corrected",(unsigned long long)corrected,"total",(unsigned long long)total)));
        Object timing(checked(Py_BuildValue("{s:d,s:d,s:d,s:d,s:d,s:d}","read_seconds",read_seconds,"shard_read_cpu_seconds",read_cpu,"filter_seconds",filter_seconds,"filter_cpu_seconds",filter_cpu,"assign_seconds",seconds(assigned),"assign_cpu_seconds",cpu()-ca)));
        return PyTuple_Pack(6,a.value,b.value,c.value,d.value,audit.value,timing.value);
    }catch(const std::exception& error){if(!PyErr_Occurred())PyErr_SetString(PyExc_ValueError,error.what());return nullptr;}
}
}
