import unittest
from pyteomics.parser import *
from string import uppercase
import random
class ParserTest(unittest.TestCase):
    def setUp(self):
        self.simple_sequences = [''.join([random.choice(uppercase) for i in range(
            int(random.uniform(1, 20)))]) for j in range(10)]

    def test_parse_sequence_simple(self):
        for seq in self.simple_sequences:
            self.assertEqual(seq, ''.join(parse_sequence(seq, labels=uppercase)))

if __name__ == '__main__':
    unittest.main()
