#!/usr/bin/env python3
"""Bounded public-source audit, no publication code executed and no large dataset downloads."""
import io,json,struct,re,urllib.request,urllib.error,hashlib,html.parser,collections,zipfile
from pathlib import Path
NATURE="https://www.nature.com/articles/s41928-026-01642-z"
IVY="https://zenodo.org/records/15095523/files/Plant_data.zip?download=1"
Z3="https://zenodo.org/records/22081982"
headers={"User-Agent":"Mozilla/5.0 (compatible; 901-source-audit/1.0)","Accept":"*/*"}
report={"stage":"E117","project":"901","scope":"bounded remote source inventory, NOT actual Bombadil physiologic forecast","targets":[]}
def request(url,method="GET",r=None,limit=800000,timeout=30):
    head=dict(headers)
    if r is not None:head["Range"]=r
    req=urllib.request.Request(url,headers=head,method=method)
    try:
        with urllib.request.urlopen(req,timeout=timeout) as resp:
            st=resp.status
            c=dict(resp.headers.items())
            body=b"" if method=="HEAD" else resp.read(limit+1)
            final=resp.geturl()
        if len(body)>limit:
            return {"status":st,"final_host":urllib.parse.urlsplit(final).hostname,"headers":c,"limit_reached":True},b""
        return {"status":st,"final_host":urllib.parse.urlsplit(final).hostname,"headers":c,"bytes_read":len(body)},body
    except Exception as ex:
        return {"error":type(ex).__name__,"reason":str(ex)[:240]},b""
class Parser(html.parser.HTMLParser):
    def __init__(self):
        super().__init__()
        self.links=[];self.current=None
    def handle_starttag(self,tag,attrs):
        if tag=="a":
            self.current={"href":dict(attrs).get("href",""),"text":""}
        elif tag in ("img","span") and self.current is not None:
            self.current["text"]+=" "
    def handle_data(self,txt):
        if self.current is not None:self.current["text"]+=txt
    def handle_endtag(self,tag):
        if tag=="a" and self.current is not None:
            self.links.append(self.current);self.current=None
t,raw=request(NATURE,limit=3500000)
nature={"site":NATURE,"fetch":{"status":t.get("status"),"error":t.get("error"),"reason":t.get("reason"),"bytes":len(raw),"content_type":t.get("headers",{}).get("Content-Type")}}
if raw and (b"nature.com" in raw[:100000] or b"Source Data" in raw):
    p=Parser()
    try:p.feed(raw.decode("utf-8","replace"))
    except Exception as ex:nature["parse_error"]=str(ex)
    links=p.links
    subset=[l for l in links if "source" in l["text"].lower() or "fig. 8" in l["text"].lower() or ".xlsx" in l["href"].lower() or "MOESM" in l["href"]]
    nature["links_relevant"]=[{"text":x["text"][:115],"href":x["href"][:290]} for x in subset[:80]]
    nature["all_href_count"]=len(links)
else:nature["links_relevant"]=[]
report["targets"].append(nature)
# Ivy 443MB: do not fetch entire archive; verify HTTP partial support, parse ZIP central metadata
head,_=request(IVY,method="HEAD",limit=0)
ivy={"site":IVY,"published_md5":"da4711e3ee553c172d13977f7ea505d9","head":{"status":head.get("status"),"reason":head.get("reason"),"content_length":head.get("headers",{}).get("Content-Length"),"accept_ranges":head.get("headers",{}).get("Accept-Ranges")}}
tail_meta,tail=request(IVY,r="bytes=-131072",limit=140000)
ivy["range_probe"]={"status":tail_meta.get("status"),"reason":tail_meta.get("reason"),"content_range":tail_meta.get("headers",{}).get("Content-Range"),"bytes":len(tail),"aborted":tail_meta.get("limit_reached",False)}
if tail_meta.get("status")==206:
    end=tail.rfind(b"PK\x05\x06")
    if end>=0 and len(tail)-end>=22:
        v=struct.unpack_from("<IHHHHIIH",tail,end)
        cd_count=v[4];cd_bytes=v[5];cd_offset=v[6]
        ivy["EOCD"]={"total_entries":cd_count,"directory_bytes":cd_bytes,"directory_offset":cd_offset}
        if cd_bytes<1200000 and cd_offset!=0xffffffff:
            cmeta,cdata=request(IVY,r=f"bytes={cd_offset}-{cd_offset+cd_bytes-1}",limit=1200000)
            ivy["central_fetch"]={"status":cmeta.get("status"),"bytes":len(cdata),"content_range":cmeta.get("headers",{}).get("Content-Range")}
            if cmeta.get("status")==206 and len(cdata)==cd_bytes:
                off=0;files=[]
                while off+46<=len(cdata) and cdata[off:off+4]==b"PK\x01\x02":
                    hdr=struct.unpack_from("<IHHHHHHIIIHHHHHII",cdata,off)
                    compression=hdr[4];size_c=hdr[8];size_u=hdr[9];nl=hdr[10];xl=hdr[11];cl=hdr[12]
                    name=cdata[off+46:off+46+nl].decode("utf-8","replace")
                    files.append({"name":name[:250],"compressed_bytes":size_c,"uncompressed_bytes":size_u,"compression":compression,"local_offset":hdr[16]})
                    off+=46+nl+xl+cl
                ivy["central_listing"]={"parsed_entries":len(files),"trailing_bytes":len(cdata)-off,"first_members":files[:50],
                    "total_compressed_listed":sum(x["compressed_bytes"] for x in files),
                    "name_extensions":dict(collections.Counter(x["name"].split(".")[-1].lower() for x in files)),
                    "filename_mentions_of_leaf_water_potential_or_biomarker":[x["name"] for x in files if re.search(r"water[_ -]?potential|biomarker|turgor|stomata|phy|fluorescence|wet|label",x["name"],re.I)][:60]}
                # For any plain text metadata <0.3MB, retrieve EXACT compressed member range only
                small=[x for x in files if x["compressed_bytes"]<180000 and x["uncompressed_bytes"]<400000 and re.search(r"readme|metadata|info|manifest|data[_-]?dictionary",x["name"],re.I)]
                ivy["small_meta_candidates"]=small[:20]
                for x in small[:2]:
                    local_off=x["local_offset"]
                    h,dat=request(IVY,r=f"bytes={local_off}-{local_off+30+x['compressed_bytes']+1000}",limit=200000)
                    if h.get("status")!=206 or dat[:4]!=b"PK\x03\x04":continue
                    lo=struct.unpack_from("<IHHHHHIIIHH",dat,0);npos=30+lo[9]+lo[10]
                    compressed=dat[npos:npos+x["compressed_bytes"]]
                    try:
                        if x["compression"]==0:payload=compressed
                        elif x["compression"]==8:
                            import zlib
                            payload=zlib.decompress(compressed,-15)
                        else:continue
                        ivy.setdefault("inspected_small_members",[]).append({"name":x["name"],"size":len(payload),"sha256":hashlib.sha256(payload).hexdigest(),"preview":payload[:800].decode("utf-8","replace")})
                    except Exception as ex:ivy.setdefault("small_member_errors",[]).append({"name":x["name"],"error":str(ex)[:120]})
report["targets"].append(ivy)
# direct 3rd-party source public GitHub scope is evaluated separately by E117_fig4_source_audit
target=Path("901_E117_publication_and_ivy_source_gate.json")
target.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf8")
for it in report["targets"]:
    if "EOCD" in it:
        brief={"EOCD":it["EOCD"],"central":it.get("central_listing",{}),"range":it["range_probe"]}
        brief["central"]= {k:v for k,v in brief["central"].items() if k!="first_members"}
        print("E117_IVY",json.dumps(brief,ensure_ascii=False)[:15000],flush=True)
    else:print("E117_NATURE",json.dumps(it,ensure_ascii=False)[:17000],flush=True)
print("E117_GATE_SHA256",hashlib.sha256(target.read_bytes()).hexdigest(),flush=True)
