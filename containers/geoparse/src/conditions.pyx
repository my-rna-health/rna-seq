import GEOparse
import pandas as pd
import pprint
import requests
from functional import seq

def extract(id: str):
    print("extraction id" % id)
    gse = GEOparse.get_GEO(id)
    pprint.pprint(vars(gse))

def download_file(url: str):
    local_filename = url.split('/')[-1]
    # NOTE the stream=True parameter
    r = requests.get(url, stream=True)
    with open(local_filename, 'wb') as f:
        for chunk in r.iter_content(chunk_size=1024):
            if chunk: # filter out keep-alive new chunks
                f.write(chunk)
                #f.flush() commented by recommendation from J.F.Sebastian
    return local_filename
