import gzip
import importlib.util
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('preflight', Path(__file__).resolve().parents[1]/'cluster/dna_preflight.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class FastqChecks(unittest.TestCase):
    def test_plain_and_compressed_prefix(self):
        with tempfile.TemporaryDirectory() as directory:
            for suffix in ['.fastq', '.fastq.gz']:
                path = Path(directory)/('reads'+suffix)
                opener = gzip.open if suffix.endswith('.gz') else open
                with opener(path, 'wt') as stream:
                    stream.write('@r1\nACGT\n+\nIIII\n@r2\nAC\n+\nII\n')
                self.assertEqual(module.fastq_prefix(path)['length_distribution'], {4: 1, 2: 1})
                self.assertEqual(module.fastq_prefix(path, limit=1)['records_examined'], 1)

    def test_truncated_fastq_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'bad.fastq'
            path.write_text('@r\nACGT\n+\n')
            with self.assertRaises(ValueError):
                module.fastq_prefix(path)


if __name__ == '__main__':
    unittest.main()
