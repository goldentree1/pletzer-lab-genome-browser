import { buildConfig } from './build-config';

const conf = buildConfig({
  ncbiName: 'GCF_031932345.1',
  firstRegion: 'NZ_JAVSCP010000001.1',
  data: {
    refSeq: 'refseq.fna.gz',
    genomic: 'genomic.gff.sorted.noregion.gff.gz',
    coverage: [],
  },
});

export default conf;
