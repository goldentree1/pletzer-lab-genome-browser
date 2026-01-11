import { gen } from './config-gen';

const g = gen({
  ncbiName: 'GCF_000013465.1',
  firstRegion: 'NC_007793.1',
  data: {
    refSeq: 'refseq.fna.gz',
    genomic: 'genomic.gff.sorted.noregion.gff.gz',
    coverage: [],
  },
});

export default g;
