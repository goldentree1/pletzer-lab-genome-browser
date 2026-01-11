import { gen } from './config-gen';

const g = gen({
  ncbiName: 'GCF_000013425.1',
  firstRegion: 'NC_007795.1',
  data: {
    refSeq: 'refseq.fna.gz',
    genomic: 'genomic.gff.sorted.noregion.gz',
    coverage: [],
  },
});

export default g;
