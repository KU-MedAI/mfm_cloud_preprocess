"""오류 응답 혼입·구간 불일치·버전 변경 및 실제 BAM/CRAM 손상 복구를 검증한다."""
import json
import random
import subprocess
import tempfile
import threading
import unittest
from http.server import BaseHTTPRequestHandler,ThreadingHTTPServer
from pathlib import Path
from unittest.mock import patch
import pysam
from handoff import http_download as http
from handoff import alignment_integrity as integrity
from handoff.alignment_recovery import repair

XML=b'<Error><Code>ConnectionClosedException</Code><Message>Premature end of Content-Length delimited message body</Message></Error>'

class LocalRemote:
    def __init__(self,data,mode='normal'):
        self.data=data;self.mode=mode;self.requests=[];self.failed_once=False
        owner=self
        class Handler(BaseHTTPRequestHandler):
            def log_message(self,*args):pass
            def do_GET(self):
                start,end=map(int,self.headers['Range'].removeprefix('bytes=').split('-'))
                owner.requests.append((start,end))
                if owner.mode=='404':self.send_error(404);return
                data=owner.data[start:end+1];range_start=start;etag='"v1"'
                if (start,end)!=(0,0) and not owner.failed_once:
                    owner.failed_once=True
                    if owner.mode=='xml':data=data[:40]+XML+data[40+len(XML):]
                    if owner.mode=='range':range_start=start+1
                    if owner.mode=='etag':etag='"v2"'
                self.send_response(206);self.send_header('ETag',etag)
                self.send_header('Last-Modified','Mon, 01 Jan 2024 00:00:00 GMT')
                self.send_header('Content-Range',f'bytes {range_start}-{end}/{len(owner.data)}')
                self.send_header('Content-Length',str(len(data)));self.end_headers()
                if owner.mode=='short' and (start,end)!=(0,0) and len(owner.requests)==2:
                    self.wfile.write(data[:10]);self.close_connection=True
                else:self.wfile.write(data)
        self.server=ThreadingHTTPServer(('127.0.0.1',0),Handler)
        self.thread=threading.Thread(target=self.server.serve_forever,daemon=True)
        self.url=f'http://127.0.0.1:{self.server.server_port}/input'
    def __enter__(self):self.thread.start();return self
    def __exit__(self,*args):self.server.shutdown();self.server.server_close();self.thread.join()

class RangeTests(unittest.TestCase):
    def test_빈파일이나본문없는과거버전표지는새전송으로복구(self):
        for exists in [False, True]:
            with self.subTest(empty_file_exists=exists), tempfile.TemporaryDirectory() as d,LocalRemote(b'x'*1024) as remote:
                partial=Path(d)/'part'
                if exists:partial.write_bytes(b'')
                metadata=Path(str(partial)+'.http.json')
                metadata.write_text(json.dumps({'url':remote.url,'size':1024,'etag':'"old"','modified':'old'}))
                http.transfer(remote.url,partial,lambda _:False,chunk_bytes=512)
                self.assertEqual(partial.read_bytes(),remote.data)
                self.assertEqual(json.loads(metadata.read_text())['etag'],'"v1"')

    def test_본문이있으면원격버전불일치를계속거부하고보존(self):
        with tempfile.TemporaryDirectory() as d,LocalRemote(b'x'*1024) as remote:
            partial=Path(d)/'part';partial.write_bytes(b'x'*256)
            metadata=Path(str(partial)+'.http.json')
            old={'url':remote.url,'size':1024,'etag':'"old"','modified':'old'}
            metadata.write_text(json.dumps(old))
            with self.assertRaises(http.RemoteChanged):http.transfer(remote.url,partial,lambda _:False)
            self.assertEqual(partial.read_bytes(),b'x'*256)
            self.assertEqual(json.loads(metadata.read_text()),old)
            self.assertEqual(remote.requests,[(0,0)])

    def test_xml_length_correct_success_response_is_retried_before_append(self):
        data=bytes(range(256))*4
        with tempfile.TemporaryDirectory() as d,LocalRemote(data,'xml') as remote,patch.object(http.time,'sleep'):
            p=Path(d)/'input.part'
            http.transfer(remote.url,p,lambda _:False,chunk_bytes=512)
            self.assertEqual(p.read_bytes(),data)
            self.assertEqual(remote.requests.count((0,511)),2)

    def test_wrong_content_range_is_not_saved(self):
        data=b'x'*1024
        with tempfile.TemporaryDirectory() as d,LocalRemote(data,'range') as remote,patch.object(http.time,'sleep'):
            p=Path(d)/'part';http.transfer(remote.url,p,lambda _:False,chunk_bytes=512)
            self.assertEqual(p.read_bytes(),data)
            self.assertEqual(remote.requests.count((0,511)),2)

    def test_incomplete_body_retries_same_range(self):
        data=b'x'*1024
        with tempfile.TemporaryDirectory() as d,LocalRemote(data,'short') as remote,patch.object(http.time,'sleep'):
            p=Path(d)/'part';http.transfer(remote.url,p,lambda _:False,chunk_bytes=512)
            self.assertEqual(p.read_bytes(),data)
            self.assertEqual(remote.requests.count((0,511)),2)

    def test_version_change_stops_before_appending(self):
        with tempfile.TemporaryDirectory() as d,LocalRemote(b'x'*1024,'etag') as remote:
            p=Path(d)/'part'
            with self.assertRaises(http.RemoteChanged):http.transfer(remote.url,p,lambda _:False)
            self.assertFalse(p.exists())

    def test_existing_partial_with_xml_is_preserved(self):
        data=b'x'*2048
        with tempfile.TemporaryDirectory() as d,LocalRemote(data) as remote:
            p=Path(d)/'part';bad=b'x'*100+XML+b'x'*100;p.write_bytes(bad)
            with self.assertRaises(http.InvalidResponse):http.transfer(remote.url,p,lambda _:False)
            self.assertEqual(p.read_bytes(),bad)

    def test_scan_detects_xml_across_read_boundaries(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'part';p.write_bytes(b'x'*255+XML+b'x'*50)
            self.assertEqual(http.error_offsets(p,chunk_bytes=256),[255])

    def test_existing_good_partial_checks_and_resumes(self):
        data=bytes(range(256))*8
        with tempfile.TemporaryDirectory() as d,LocalRemote(data) as remote:
            p=Path(d)/'part';p.write_bytes(data[:1000]);http.transfer(remote.url,p,lambda _:False,chunk_bytes=512)
            self.assertEqual(p.read_bytes(),data)
            self.assertIn((1000,1511),remote.requests)

class AlignmentTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.bam=self.root/'good.bam';self.ref=self.root/'ref.fa'
        self.ref.write_text('>chr1\n'+'A'*10000+'\n');pysam.faidx(str(self.ref))
        rand=random.Random(71)
        header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr1','LN':10000}],'RG':[{'ID':'r','SM':'S'}]}
        with pysam.AlignmentFile(str(self.bam),'wb',header=header) as stream:
            for i in range(1500):
                r=pysam.AlignedSegment();r.query_name='read_'+str(i)+'_'+str(rand.getrandbits(128))
                r.query_sequence=''.join(rand.choices('ACGT',k=100));r.query_qualities=pysam.qualitystring_to_array('I'*100)
                r.flag=0;r.reference_id=0;r.reference_start=i;r.mapping_quality=60;r.cigarstring='100M';r.set_tag('RG','r')
                stream.write(r)
    def tearDown(self):self.tmp.cleanup()

    def test_bam_all_blocks_and_fingerprint_cache(self):
        report=self.root/'physical.json'
        result=integrity.validate(self.bam,report)
        self.assertEqual(result['status'],'pass');self.assertGreater(result['blocks_checked'],1)
        raw=bytearray(self.bam.read_bytes());raw[1000]^=0xff;self.bam.write_bytes(raw)
        with self.assertRaises(integrity.AlignmentCorruption):integrity.validate(self.bam,report)
        self.assertEqual(json.loads(report.read_text())['status'],'failed')

    def test_actual_cram_all_container_and_block_crc(self):
        for version in ['3.0','3.1']:
            with self.subTest(version=version):
                cram=self.root/('good_'+version+'.cram')
                with pysam.AlignmentFile(str(self.bam),'rb') as source,pysam.AlignmentFile(str(cram),'wc',header=source.header,reference_filename=str(self.ref),format_options=[('version='+version).encode()]) as target:
                    for r in source:target.write(r)
                result=integrity.scan(cram)
                self.assertEqual(result['format'],'CRAM'+version);self.assertGreater(result['blocks_checked'],1)
                raw=bytearray(cram.read_bytes());raw[-45]^=0xff;cram.write_bytes(raw)
                with self.assertRaises(integrity.AlignmentCorruption):integrity.scan(cram)

    def test_wgs_replacement_updates_resource_index_and_bound_manifest(self):
        from handoff import resume,cli
        row={'sample_id':'HG00151','alignment_path':'wgs/HG00151.cram','rna_path':'rna/HG00151.bam',
             'platform':'illumina','project_id':'1KGP_ILLUMINA','public_vcf_path':'','data_group':'wgs_rna_pair'}
        config={'data_root':str(self.root),'state_root':str(self.root/'state'),'output_root':str(self.root/'out'),
                'run_id':'test','sample_wgs_paths':{'HG00151':'recovery/HG00151.cram'}}
        resources=[{'kind':'wgs_alignment','path':row['alignment_path'],'url':'https://example.test/original.cram'},
                   {'kind':'wgs_index','path':row['alignment_path']+'.crai','url':'https://example.test/original.cram.crai'}]
        result=resume.sample_downloads(config,row,resources)
        self.assertEqual(result[0]['path'],'recovery/HG00151.cram')
        self.assertEqual(result[1]['path'],'recovery/HG00151.cram.crai')
        self.assertEqual(resources[0]['path'],'wgs/HG00151.cram')
        self.assertEqual(result[0]['url'],resources[0]['url'])
        refs={'ncbi':'ncbi.fa','autosomes':'auto.fa','broad':'broad.fa','gtf':'genes.gtf','rna_extra':'extra.fa'}
        with patch.object(cli.planner,'read_tsv',return_value=[]),patch.object(cli.bundled_assets,'panel_directory',return_value=self.root/'vcf'),patch('handoff.references.exact_cram_reference',return_value=self.ref):
            _,bound=cli.bind_manifest(config,row,result,refs)
        self.assertEqual(bound['alignment_path'],str(self.root/'recovery/HG00151.cram'))
        config['sample_wgs_paths']['HG00151']='../outside.cram'
        with self.assertRaises(ValueError):resume.sample_downloads(config,row,resources)

    def test_http404_preserves_index_fallback_exception_contract(self):
        from handoff import download_retry
        with LocalRemote(b'x'*1024,'404') as remote:
            with self.assertRaises(subprocess.CalledProcessError) as error:
                download_retry.transfer(remote.url,self.root/'missing.part',lambda _:False)
            self.assertEqual(error.exception.returncode,22)
            self.assertEqual(len(remote.requests),1)

    def test_repair_preserves_bad_source_and_restores_full_bam(self):
        good=self.bam.read_bytes();bad=bytearray(good);bad[1000:1000+len(XML)]=XML
        self.bam.write_bytes(bad);before=self.bam.read_bytes();target=self.root/'repaired.bam'
        with LocalRemote(good) as remote:
            result=repair(self.bam,target,remote.url,self.root/'recovery.json',len(good),'"v1"')
        self.assertEqual(result['status'],'pass');self.assertEqual(target.read_bytes(),good)
        self.assertEqual(self.bam.read_bytes(),before)
        self.assertEqual(integrity.scan(target)['bytes_checked'],len(good))

    def test_missing_bam_eof_is_rejected(self):
        self.bam.write_bytes(self.bam.read_bytes()[:-28])
        with self.assertRaises(integrity.AlignmentCorruption):integrity.scan(self.bam)

if __name__=='__main__':unittest.main()
