import { gen } from './config-gen';

const g = gen({
  ncbiName: 'GCF_000006765.1',
  firstRegion: 'NC_002516.2',
  data: {
    refSeq: 'refseq.fna.gz',
    genomic: 'genomic.gff.sorted.noregion.gff.gz',
    coverage: [],
  },
});

export default g;
