import { gen } from './config-gen';

const genConf = gen({
  ncbiName: 'GCF_000014625',
  firstRegion: 'NC_008463.1',
  data: {
    refSeq: 'GCF_000014625.1_ASM1462v1_genomic.fna.gz',
    genomic: 'genomic.gff.sorted.noregion.gff.gz',
    coverage: ['PA14_I_1.bw', 'PA14_Un_1.bw'],
  },
});

export default genConf;
