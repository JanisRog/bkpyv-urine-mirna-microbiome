import importlib.util
from pathlib import Path
import tempfile
import unittest

spec=importlib.util.spec_from_file_location('recount',Path(__file__).resolve().parents[1]/'cluster/recount_mirna.py')
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)


def arf(query, precursor):
    return '\t'.join([query,'20','1','20','a'*20,precursor,'20','10','29','a'*20,'+','0','m'*20])+'\n'


class RecountTests(unittest.TestCase):
    def test_same_mature_and_distinct_mature_ambiguity(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp)
            e=d/'expression.tsv'; mature=d/'mature.arf'; reads=d/'reads.arf'
            e.write_text('#miRNA\tread_count\tprecursor\nA\t8\tp1\nA\t8\tp2\nB\t3\tp3\n')
            mature.write_text(arf('A','p1')+arf('A','p2')+arf('B','p3'))
            reads.write_text(arf('hsa_1_x5','p1')+arf('hsa_1_x5','p2')+
                             arf('hsa_2_x3','p1')+arf('hsa_2_x3','p2')+arf('hsa_2_x3','p3'))
            _,h,u,unique,f,groups,q=m.recount(e,mature,reads)
            self.assertEqual(dict(h),{'A':16,'B':3})
            self.assertEqual(dict(u),{'A':8,'B':3})
            self.assertEqual(dict(unique),{'A':5})
            self.assertEqual(dict(f),{'A':6.5,'B':1.5})
            self.assertEqual(groups[('A','B')],3)
            self.assertEqual(q['distinct_assigned_read_total'],8)
            e.write_text(e.read_text().replace('B\t3','B\t4'))
            with self.assertRaisesRegex(ValueError,'reconstruction failed'):
                m.recount(e,mature,reads)


if __name__=='__main__': unittest.main()
