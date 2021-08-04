#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Extract methylation from fast5 files into a RocksDB file. Also has an interface to read those values.

Created on Thursday, 25. July 2019.
"""
import glob
import os.path


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Extract methylation from fast5 files")

    parser.add_argument("input_fast5_dir", type=str,
                        help="Input dir of fast5 files  [default:%(default)s]")
    parser.add_argument("-d", "--mod_data",
                        help="Database to store the modifications to  [default:%(default)s]",
                        default="base_mods.rocksdb")
    parser.add_argument("-p", "--processes", type=int,
                        help="Database to store the modifications to  [default:%(default)s]",
                        default=1)
    parser.add_argument("-V", "--verbose", default=False, action="store_true",
                        help="Be more verbose with output [default:%(default)s]")

    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s:%(funcName)s:%(levelname)s:%(message)s')

    return args


class MethylDB(object):
    def __init__(self, db_name, q=None):
        import rocksdb
        self._db_name = db_name
        self._q = q

        if self._q is None:
            opts = rocksdb.Options()
            opts.create_if_missing = True
            opts.max_open_files = 300000
            opts.max_open_files = -1  # Dangerous
            opts.write_buffer_size = 2 * 512 * 1024 ** 2
            opts.max_write_buffer_number = 3
            opts.target_file_size_base = 512 * 1024 ** 2  # MB

            # opts.compression = rocksdb.CompressionType.zlib_compression

            opts.table_factory = rocksdb.BlockBasedTableFactory(
                    filter_policy=rocksdb.BloomFilterPolicy(10),
                    # cache_index_and_filter_blocks=True,
                    # optimize_filters_for_hits=True,
                    block_cache=rocksdb.LRUCache(5 * (1024 ** 3)),
                    block_size=64 * 1024,
                    block_cache_compressed=rocksdb.LRUCache(500 * (1024 ** 2)))

            self._db = rocksdb.DB(self._db_name, opts)
        A, mA, C, mC, G, T = 0, 1, 2, 3, 4, 5

    def close(self):
        del (self._db)

    def _put(self, *args):
        if self._q:
            self._q.put(args)
        else:
            self._db.put(*args)

    def __len__(self):
        "Return approximate number of key-value pairs in the database"
        return int(self._db.get_property(b"rocksdb.estimate-num-keys"))

    def put(self, read_id, likelihoods, sequence=None):
        from uuid import UUID
        from numpy import uint8, ndarray

        assert isinstance(likelihoods, ndarray), "Likelihoods must be of type ndarray"
        assert likelihoods.ndim == 1, "Likelihoods must be one dimensional"
        assert likelihoods.dtype == uint8, "Likelihoods must be of dtype uint8"

        read_uuid = UUID(read_id)

        self._put(read_uuid.bytes, likelihoods.tobytes())
        # else:
        #    self._db.put(*args)

        if sequence is not None:
            self._put(read_uuid.bytes + b"/seq", sequence.encode("ascii"))

    def get(self, read_id, with_sequence=False):
        import numpy as np
        from uuid import UUID

        read_uuid = UUID(read_id)
        mod_data = self._db.get(read_uuid.bytes)
        if mod_data is None:
            likelihoods = None
        else:
            likelihoods = np.frombuffer(mod_data, dtype=np.uint8)

        if with_sequence:
            likelihoods = likelihoods, self._db.get(read_uuid.bytes + b"/seq")
        return likelihoods

    def update_fast5(self, fast5_filepath, mod_index=3, verbose=False):
        """Update (i.e. add or change) the methylation data for reads in the given fast5 file.
        
        mod_index gives the index of the modification call table to store in the database. 
                    Default is mC modification. Indices: A,mA,C,mC,G,T = 0,1,2,3,4,5"""

        from ont_fast5_api.fast5_interface import get_fast5_file
        import numpy as np
        import logging as log
        if verbose:
            from tqdm import tqdm
        else:
            def tqdm(x):
                return x

        log.info("Processing file {}".format(fast5_filepath))

        UNMODIFIED_BASES = [b"A", b"A", b"C", b"C", b"G", b"T"]
        assert mod_index >= 0 and mod_index < len(UNMODIFIED_BASES), "mod_index must be in the range 0-5."

        BASE = UNMODIFIED_BASES[mod_index]

        log.info("Looking for modification {} of base {}.".format(mod_index, BASE))

        with get_fast5_file(fast5_filepath, mode="r") as f5:
            for read_id in tqdm(f5.get_read_ids()):
                # if read_idx%100:
                #    log.info("Processing read {}".format(read_id))

                read = f5.get_read(read_id)

                latest_basecall = read.get_latest_analysis('Basecall_1D')

                mod_base_table = read.get_analysis_dataset(
                        latest_basecall, 'BaseCalled_template/ModBaseProbs')
                if mod_base_table is None:
                    log.info("No ModBaseProbs for {}".format(read_id))
                    continue

                fastq = read.get_analysis_dataset(
                        latest_basecall, 'BaseCalled_template/Fastq')
                if fastq is None:
                    log.info("No Fastq for {}".format(read_id))
                    continue

                seq_title, seq, _, qvals, _ = fastq.split("\n")

                mod_likelihoods = mod_base_table[np.fromstring(seq, "|S1") == BASE, mod_index]

                self.put(read_id, mod_likelihoods)

                # assert (self.get(read_id) == mod_likelihoods).all(),"Mismatch on "+read_id


def _fast5_putter(fname, q):
    import logging as log
    import os
    log.info("Processing {} to {} in process {}.".format(fname, q, os.getpid()))
    mdb = MethylDB(None, q=q)
    mdb.update_fast5(fname)
    return (fname, str(q))


if __name__ == '__main__':
    args = main()
    mdb = MethylDB(args.mod_data)

    import logging as log

    log.info(args)

    indir = args.input_fast5_dir
    fnlist = glob.glob(os.path.join(indir, '*.fast5'))
    log.info(f'Total files={len(fnlist)}')

    if args.processes == 1:
        for fn in fnlist:
            mdb.update_fast5(fn, verbose=args.verbose)
    elif args.processes > 1:
        from tqdm import tqdm
        import multiprocessing as mp
        import itertools as it
        from queue import Empty

        procs = min(args.processes, len(fnlist))
        log.info("Will run {} processes in parallel.".format(procs))

        m = mp.Manager()
        q = m.Queue(1000)

        with mp.Pool(procs) as pool:
            # Read the fast5 files in parallel and put the results in queue q
            read_result = pool.starmap_async(_fast5_putter, zip(fnlist, it.repeat(q)))
            for _ in tqdm(it.repeat(True)):
                try:
                    data = q.get(timeout=1)
                except Empty:
                    if read_result.ready():
                        # log.info("Fast5 processing is ready. Got {}".format(read_result.get()))
                        log.info("Fast5 processing is ready.")
                        break
                    else:
                        log.info("Stalling... Fast5 processing takes time.")
                else:
                    mdb._db.put(*data)
