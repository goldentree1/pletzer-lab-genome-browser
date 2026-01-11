import { buildConfig } from './build-config';

const conf = buildConfig({
  ncbiName: 'GCF_000281535.2',
  firstRegion: 'NZ_CP008827.1',
  data: {
    refSeq: 'refseq.fna.gz',
    genomic: 'genomic.gff.sorted.noregion.gff.gz',
    coverage: [],
  },
});

export default conf;
