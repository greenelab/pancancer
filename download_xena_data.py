# This script will download publicly available pancancer data that
# does not have restricted use
# Modified from https://github.com/cognoma/machine-learning/1.download.ipynb

import os
from urllib import request
import requests
from subprocess import call

if not os.path.exists('data/xena'):
    os.makedirs('data/xena')

figshare_id = "3487685"

# Use Figshare API to download UCSC xena data
fig_url = "https://api.figshare.com/v2/articles/{}".format(figshare_id)
response = requests.get(fig_url)
response_json = response.json()

for file_data in response_json['files']:
    file_url = file_data['download_url']
    file_name = file_data['name']
    path = os.path.join('data', 'xena', file_name)
    if not os.path.exists(path) and file_name != 'mutation-matrix.tsv.bz2':
        request.urlretrieve(file_url, path)

# Download gene maps from cognoma/cancer-data
gene_map = 'https://raw.githubusercontent.com/cognoma/cancer-data/'\
           '54140cf6addc48260c9723213c40b628d7c861da/mapping/HiSeqV2-genes/'\
           'HiSeqV2-gene-map.tsv'

if not os.path.exists('data/xena/HiSeqV2-gene-map.tsv'):
    request.urlretrieve(gene_map, os.path.join('data', 'xena',
                                               'HiSeqV2-gene-map.tsv'))

original_mutation = 'https://zenodo.org/record/56735/files/'\
                    'gbm_classifier_data.tar.gz'

if not os.path.exists('data/xena/gbm_classifier_data.tar.gz'):
    request.urlretrieve(original_mutation,
                        os.path.join('data', 'xena',
                                     'gbm_classifier_data.tar.gz'))

    call(['tar', '-zxvf', 'data/xena/gbm_classifier_data.tar.gz', '-C',
          'data/xena/'])
