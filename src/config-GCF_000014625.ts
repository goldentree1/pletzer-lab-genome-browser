import type { JBrowseConfig } from './types';

const config: JBrowseConfig = {
  assembly: {
    name: 'GCF_000014625',
    sequence: {
      type: 'ReferenceSequenceTrack',
      trackId: 'RefSeqTrack',
      adapter: {
        type: 'BgzipFastaAdapter',
        uri: '/data/GCF_000014625/GCF_000014625.1_ASM1462v1_genomic.fna.gz',
      },
    },
  },
  tracks: [
    {
      type: 'FeatureTrack',
      trackId: 'GFF3GeneTrack',
      name: 'GFF3 Track',
      assemblyNames: ['GCF_000014625'],
      adapter: {
        type: 'Gff3TabixAdapter',
        gffGzLocation: {
          uri: '/data/GCF_000014625/genomic.gff.sorted.noregion.gff.gz',
          locationType: 'UriLocation',
        },
        index: {
          location: {
            uri: '/data/GCF_000014625/genomic.gff.sorted.noregion.gff.gz.tbi',
            locationType: 'UriLocation',
          },
        },
      },
    },
  ],
  defaultSession: {
    name: 'this session',
    margin: 0,
    view: {
      id: 'linearGenomeView',
      type: 'LinearGenomeView',
      displayedRegions: [
        {
          refName: 'NC_008463.1',
          start: 0,
          end: 6537648,
          assemblyName: 'GCF_000014625',
        },
      ],
      init: {
        assembly: 'GCF_000014625',
        loc: 'NC_008463.1:1..5,000',
        tracks: ['RefSeqTrack', 'GFF3GeneTrack'],
      },
    },
  },

  aggregateTextSearchAdapters: [
    {
      type: 'TrixTextSearchAdapter',
      textSearchAdapterId: 'lesb58-index',
      ixFilePath: {
        uri: '/data/GCF_000014625/trix/GCF_000014625.ix',
        locationType: 'UriLocation',
      },
      ixxFilePath: {
        uri: '/data/GCF_000014625/trix/GCF_000014625.ixx',
        locationType: 'UriLocation',
      },
      metaFilePath: {
        uri: '/data/GCF_000014625/trix/GCF_000014625_meta.json',
        locationType: 'UriLocation',
      },
      assemblyNames: ['GCF_000014625'],
    },
  ],
};

export default config;
