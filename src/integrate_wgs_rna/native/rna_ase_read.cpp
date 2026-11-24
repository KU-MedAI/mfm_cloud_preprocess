// QNAME별 입력 순서와 동점의 첫 mate를 보존하며 shard를 Python 객체로 직접 해석함.
#include <Python.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <cstdint>
#include <cstring>
#include <stdexcept>

struct Object {
    PyObject* value;
    explicit Object(PyObject* v=nullptr):value(v){}
    ~Object(){Py_XDECREF(value);}
    Object(const Object&)=delete;Object& operator=(const Object&)=delete;
    PyObject* release(){PyObject* v=value;value=nullptr;return v;}
};
struct Input {
    int fd=-1;size_t size=0,pos=0;const char* data=nullptr;
    ~Input(){if(data)munmap(const_cast<char*>(data),size);if(fd>=0)close(fd);}
    const char* take(size_t n){if(n>size-pos)throw std::runtime_error("잘린 ASE shard");const char* p=data+pos;pos+=n;return p;}
    uint64_t number(size_t n){const auto* p=reinterpret_cast<const unsigned char*>(take(n));uint64_t v=0;for(size_t i=0;i<n;++i)v|=uint64_t(p[i])<<(8*i);return v;}
};
static PyObject* checked(PyObject* object){if(!object)throw std::runtime_error("Python 객체 생성 실패");return object;}
static void success(int code){if(code<0)throw std::runtime_error("Python 객체 기록 실패");}
static int64_t signed64(uint64_t value){int64_t out;std::memcpy(&out,&value,sizeof(out));return out;}
static PyObject* read_shard(PyObject*,PyObject* args) {
    PyObject* rawpath;int minimum,duplicates;
    if(!PyArg_ParseTuple(args,"Oip:read_shard",&rawpath,&minimum,&duplicates))return nullptr;
    try {
        Object path;if(!PyUnicode_FSConverter(rawpath,&path.value))return nullptr;
        Input input;input.fd=open(PyBytes_AS_STRING(path.value),O_RDONLY);
        if(input.fd<0)return PyErr_SetFromErrnoWithFilenameObject(PyExc_OSError,rawpath);
        struct stat st;if(fstat(input.fd,&st)<0)return PyErr_SetFromErrno(PyExc_OSError);
        if(st.st_size<12)throw std::runtime_error("잘린 ASE shard");
        input.size=size_t(st.st_size);
        void* mapping=mmap(nullptr,input.size,PROT_READ,MAP_PRIVATE,input.fd,0);
        if(mapping==MAP_FAILED)return PyErr_SetFromErrno(PyExc_OSError);
        input.data=static_cast<const char*>(mapping);
        if(std::memcmp(input.take(8),"MFMASE1\0",8))throw std::runtime_error("ASE shard magic 오류");
        uint32_t n=uint32_t(input.number(4));if(n>1000000)throw std::runtime_error("ASE reference 개수 오류");
        Object refs(checked(PyList_New(n))),unique(checked(PyDict_New())),multi(checked(PyDict_New())),dups(checked(PySet_New(nullptr)));
        for(uint32_t i=0;i<n;++i){size_t len=size_t(input.number(4));if(len>1<<20)throw std::runtime_error("ASE reference 이름 한도 초과");PyList_SET_ITEM(refs.value,i,checked(PyUnicode_DecodeUTF8(input.take(len),Py_ssize_t(len),"strict")));}
        uint64_t duplicate_count=0;
        while(input.pos<input.size) {
            uint32_t len=uint32_t(input.number(4));uint16_t flag=uint16_t(input.number(2));int mq=int(input.number(1)),has_as=int(input.number(1));
            int64_t score=signed64(input.number(8));int32_t tid=int32_t(input.number(4));uint32_t blocks=uint32_t(input.number(4));
            if(len>1<<20||blocks>1<<20||tid<0||uint32_t(tid)>=n)throw std::runtime_error("ASE shard record 범위 오류");
            Object name(checked(PyUnicode_DecodeUTF8(input.take(len),len,"strict")));
            const auto* rawblocks=reinterpret_cast<const unsigned char*>(input.take(size_t(blocks)*16));
            if(flag&1024){++duplicate_count;success(PySet_Add(dups.value,name.value));if(!duplicates)continue;}
            if(!has_as||(mq==0&&mq<minimum))continue;
            PyObject* table=mq>=minimum?unique.value:multi.value;
            PyObject* slot=PyDict_GetItemWithError(table,name.value);if(!slot&&PyErr_Occurred())return nullptr;
            if(!slot){Object fresh(checked(PyList_New(2)));Py_INCREF(Py_None);Py_INCREF(Py_None);PyList_SET_ITEM(fresh.value,0,Py_None);PyList_SET_ITEM(fresh.value,1,Py_None);success(PyDict_SetItem(table,name.value,fresh.value));slot=PyDict_GetItem(table,name.value);}
            int mate=(flag&64)?0:1;PyObject* previous=PyList_GET_ITEM(slot,mate);
            if(previous!=Py_None&&score<=PyLong_AsLongLong(PyTuple_GET_ITEM(previous,0)))continue;
            Object values(checked(PyList_New(blocks)));
            for(uint32_t i=0;i<blocks;++i){uint64_t a=0,b=0;for(int k=0;k<8;++k){a|=uint64_t(rawblocks[i*16+k])<<(8*k);b|=uint64_t(rawblocks[i*16+8+k])<<(8*k);}PyList_SET_ITEM(values.value,i,checked(Py_BuildValue("(LL)",static_cast<long long>(signed64(a)),static_cast<long long>(signed64(b)))));}
            Object record(checked(PyTuple_New(3)));PyTuple_SET_ITEM(record.value,0,checked(PyLong_FromLongLong(score)));
            PyObject* ref=PyList_GET_ITEM(refs.value,tid);Py_INCREF(ref);PyTuple_SET_ITEM(record.value,1,ref);PyTuple_SET_ITEM(record.value,2,values.release());
            success(PyList_SetItem(slot,mate,record.release()));
        }
        Object count(checked(PyLong_FromUnsignedLongLong(duplicate_count)));
        return PyTuple_Pack(4,unique.value,multi.value,count.value,dups.value);
    }catch(const std::exception& error){if(!PyErr_Occurred())PyErr_SetString(PyExc_ValueError,error.what());return nullptr;}
}
#include "rna_ase_fused.hpp"

static PyMethodDef methods[]={{"prepare_assign",ase_fused::prepare,METH_VARARGS,"검증된 fragment 규칙의 읽기 전용 색인을 준비함."},{"assign_shards",ase_fused::run,METH_VARARGS,"두 HAP shard를 count로 직접 집계함."},{"read_shard",read_shard,METH_VARARGS,"기존 Python 해석과 같은 ASE fragment 구조를 반환함."},{nullptr,nullptr,0,nullptr}};
static PyModuleDef module={PyModuleDef_HEAD_INIT,"rna_ase_read",nullptr,-1,methods};
PyMODINIT_FUNC PyInit_rna_ase_read(){return PyModule_Create(&module);}
