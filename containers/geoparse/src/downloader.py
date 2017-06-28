import os
import re
from pprint import pprint
from typing import *

import GEOparse
from GEOparse import *
from GEOparse import utils
from functional import *


class Downloader:
    def __init__(self, folder: str = "./", temp: str = "/tmp", email: str = "antonkulaga@gmail.com"):
        self.folder = folder
        self.temp = temp
        self.email = email

    def load_series(self, geo_id: str) -> GSE:
        return GEOparse.get_GEO(geo_id, destdir=self.temp)

    def download_gse(self, gse: GSE) -> None:
        gse.download_supplementary_files(self.folder, True, email=self.email)

    def download_gsm(self, gsm_id: str, sra_filetype: str = "sra", create_folder: bool = False) -> Tuple[str, str]:
        gsm = cast(GSM, GEOparse.get_GEO(gsm_id, destdir=self.temp))
        if create_folder:
            path = os.path.join(self.folder, gsm_id)
            utils.mkdir_p(path)
            gsm.download_supplementary_files(path, True, sra_filetype=sra_filetype, email=self.email)
            return gsm_id, path
        else:
            gsm.download_supplementary_files(self.folder, True, sra_filetype=sra_filetype, email=self.email)
            directory_path = os.path.abspath(os.path.join(self.folder, "%s_%s_%s" % ('Supp',
                                                                                     gsm.get_accession(),
                                                                                     re.sub(r'[\s\*\?\(\),\.;]', '_',
                                                                                            gsm.metadata['title'][0])
                                                                                     # the directory name cannot contain many of the signs
                                                                                     )))
            return gsm_id, directory_path

    def download_samples(self, samples: List[str], sra_filetype: str = "sra", create_folder: bool = True) -> dict:
        return seq(samples).map(lambda gsm_id: self.download_gsm(gsm_id, sra_filetype, create_folder)).dict()

    # https://doc-0g-8s-docs.googleusercontent.com/docs/securesc/0s9gddvd2eanhj99krjjabvcr59io0q1/9etosuaf0i1blj7vse6q5gtaqdpnc0tk/1498060800000/17884921474303815914/17884921474303815914/0B53nZMQJitqOb0VILWRpZE5CR3c?e=download&nonce=3q9mtfvk6rcn0&user=17884921474303815914&hash=ant7ith4p3v3uvi4bf0s0smg6f0aostm
    # http://res.cloudinary.com/hrscywv4p/image/upload/c_limit,fl_lossy,h_1500,w_2000,f_auto,q_auto/v1/981253/671c739a-6648-4c02-b645-a4a8a3966cad_xnxecs.png
    def download_file(self, url: str):
        print("functono not implemented yet!")
        utils.download_from_url(url, self.folder, False )
        pass