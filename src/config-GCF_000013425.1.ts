import type { JBrowseConfig } from './types';

const config: JBrowseConfig = {
  assembly: {
    name: 'GCF_000013425.1',
    sequence: {
      type: 'ReferenceSequenceTrack',
      trackId: 'RefSeqTrack',
      adapter: {
        type: 'BgzipFastaAdapter',
        uri: '/data/GCF_000013425.1/refseq.fna.gz',
      },
    },
  },
  tracks: [
    {
      type: 'FeatureTrack',
      trackId: 'GFF3GeneTrack',
      name: 'GFF3 Track',
      assemblyNames: ['GCF_000013425.1'],
      adapter: {
        type: 'Gff3TabixAdapter',
        gffGzLocation: {
          uri: '/data/GCF_000013425.1/genomic.gff.sorted.noregion.gz',
          locationType: 'UriLocation',
        },
        index: {
          location: {
            uri: '/data/GCF_000013425.1/genomic.gff.sorted.noregion.gz.tbi',
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
          refName: 'NC_007795.1',
          start: 0,
          end: 2821361,
          assemblyName: 'GCF_000013425.1',
        },
      ],
      init: {
        assembly: 'GCF_000013425.1',
        loc: 'NC_007795.1:1..5,000',
        tracks: ['RefSeqTrack', 'GFF3GeneTrack'],
      },
    },
  },

  aggregateTextSearchAdapters: [
    {
      type: 'TrixTextSearchAdapter',
      textSearchAdapterId: 'lesb58-index',
      ixFilePath: {
        uri: '/data/GCF_000013425.1/trix/GCF_000013425.1.ix',
        locationType: 'UriLocation',
      },
      ixxFilePath: {
        uri: '/data/GCF_000013425.1/trix/GCF_000013425.1.ixx',
        locationType: 'UriLocation',
      },
      metaFilePath: {
        uri: '/data/GCF_000013425.1/trix/GCF_000013425.1_meta.json',
        locationType: 'UriLocation',
      },
      assemblyNames: ['GCF_000013425.1'],
    },
  ],
};

export default config;
