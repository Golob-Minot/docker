#!/usr/bin/env python3
import boto3
import pandas as pd
from Bio import SeqIO
import tempfile
import gzip
import io


if __name__ == "__main__":
    s3_client = boto3.client('s3')

    list_res = s3_client.list_objects_v2(
        Bucket='fh-pi-fredricks-d',
        Prefix='lab/golob/gvhd/metagenomics/20181120/qc/combined/',
        MaxKeys=1000,
    )
    combined_keys = [k['Key'] for k in list_res['Contents']]
    read_count_df = pd.DataFrame(
        columns=[
            'Bucket',
            'Key'
        ]
    )
    read_count_df['Key'] = combined_keys
    read_count_df['Bucket'] = 'fh-pi-fredricks-d'

    read_count_df['filename'] = read_count_df.Key.apply(lambda k: k.split('/')[-1])
    read_count_df['specimen'] = read_count_df.filename.apply(lambda fn: fn[15:].split('.')[0].lower())
    read_count_df['readnum'] = read_count_df.filename.apply(lambda fn: int(fn[15:].split('.')[1][1]))

    for idx, row in read_count_df.iterrows():
        print(row.specimen, row.readnum)
        with tempfile.TemporaryFile() as tfh:
            s3_client.download_fileobj(
                Fileobj=tfh,
                Bucket=row.Bucket,
                Key=row.Key,
            )
            tfh.seek(0)
            tfh_gz = io.TextIOWrapper(gzip.GzipFile(fileobj=tfh, mode='r'))
            n = 0
            for sr in SeqIO.parse(tfh_gz, 'fastq'):
                n += 1
        read_count_df.loc[idx, 'num_nonhuman_reads'] = n

    with io.BytesIO() as buf:
        with io.TextIOWrapper(buf) as text_buf:
            read_count_df.to_csv(
                text_buf,
                index=False
            )
            buf.seek(0)
            s3_client.upload_fileobj(
                Fileobj=buf,
                Bucket='fh-pi-fredricks-d',
                Key='lab/golob/gvhd/metagenomics/20181120/nonhuman_reads.csv'
            )
