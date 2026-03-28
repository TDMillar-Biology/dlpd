''' 
TDMILLAR May 2025
Parse repeat masker '.out' file
'''

class RepeatMaskerEntry:
    def __init__(self, chrom, start, end, repeat_name, class_family):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.repeat_name = repeat_name
        self.class_family = class_family

    def to_bed(self):
        return f"{self.chrom}\t{self.start}\t{self.end}\t{self.repeat_name}\t{self.class_family}"

class ParseRepeatMasker:
    @staticmethod
    def parse_repeat_masker_out(path):
        entries = []
        with open(path) as f:
            for line in f:
                if line.strip().startswith("SW") or line.strip().startswith("score") or line.strip() == "":
                    continue ## skip the headers (should be first 3 lines)
                fields = line.strip().split()
                chrom = fields[4]
                start = fields[5]
                end = fields[6]
                repeat_name = fields[9]
                class_family = fields[10]
                entries.append(RepeatMaskerEntry(chrom, start, end, repeat_name, class_family))
        return entries
