import type { JBrowseCustomConfig } from './types';

const myConfig: { [key: string]: JBrowseCustomConfig } = {
  'P.aeruginosa PAO1': {
    ncbiName: 'GCF_000006765.1',
    firstRegion: 'NC_002516.2',
    data: {
      genomic: 'genomic.gff.sorted.noregion.gff.gz',
      coverage: [],
    },
  },
  'S.aureus HG001 (NCTC 8325)': {
    ncbiName: 'GCF_000013425.1',
    firstRegion: 'NC_007795.1',
    data: {
      genomic: 'genomic.gff.sorted.noregion.gz',
      coverage: [],
    },
  },
  'S.aureus USA300LAC': {
    ncbiName: 'GCF_000013465.1',
    firstRegion: 'NC_007793.1',
    data: {
      genomic: 'genomic.gff.sorted.noregion.gff.gz',
      coverage: [],
    },
  },
  'P.aeruginosa PA14': {
    ncbiName: 'GCF_000014625',
    firstRegion: 'NC_008463.1',
    data: {
      refSeq: 'GCF_000014625.1_ASM1462v1_genomic.fna.gz',
      genomic: 'genomic.gff.sorted.noregion.gff.gz',
      coverage: [
        [
          'PA14_I_Average.bw',
          'PA14_I_1.reheadered.bw',
          'PA14_I_2.reheadered.bw',
          'PA14_I_3.reheadered.bw',
          'PA14_I_1.bin10.CPM-normalised.bw',
          'PA14_I_1.bin10.RPKM-normalised.bw',
        ],
        [
          'PA14_Un_Average.bw',
          'PA14_Un_1.reheadered.bw',
          'PA14_Un_2.reheadered.bw',
          'PA14_Un_3.reheadered.bw',
          'PA14_Un_1.bin10.CPM-normalised.bw',
          'PA14_Un_1.bin10.RPKM-normalised.bw',
        ],
      ],
    },
    extras: [],
  },
  'K.pneumoniae KPNIH1': {
    ncbiName: 'GCF_000281535.2',
    firstRegion: 'NZ_CP008827.1',
    data: {
      genomic: 'genomic.gff.sorted.noregion.gff.gz',
      coverage: [],
    },
  },
  'A.baumannii AB5075-UW': {
    ncbiName: 'GCF_031932345.1',
    firstRegion: 'NZ_JAVSCP010000001.1',
    data: {
      genomic: 'genomic.gff.sorted.noregion.gff.gz',
      coverage: [],
    },
  },
  'P.aeruginosa LESB58': {
    ncbiName: 'GCF_000026645.1',
    firstRegion: 'NC_011770.1',
    data: {
      refSeq: 'LESB58_ASM2664v1.fasta.gz',
      genomic: 'LESB58_ASM2664v1.gff.sorted.gff.noregion.gz',
      coverage: [['bigwig-bin5.bw'], ['bigwig-bin100.bw']],
    },
  },
};

export default myConfig;
