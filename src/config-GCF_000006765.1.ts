import type { JBrowseConfig } from './types';

const config: JBrowseConfig = {
  assembly: {
    name: 'GCF_000006765.1',
    sequence: {
      type: 'ReferenceSequenceTrack',
      trackId: 'RefSeqTrack',
      adapter: {
        type: 'BgzipFastaAdapter',
        uri: '/data/GCF_000006765.1/refseq.fna.gz',
      },
    },
  },
  tracks: [
    {
      type: 'FeatureTrack',
      trackId: 'GFF3GeneTrack',
      name: 'GFF3 Track',
      assemblyNames: ['GCF_000006765.1'],
      adapter: {
        type: 'Gff3TabixAdapter',
        gffGzLocation: {
          uri: '/data/GCF_000006765.1/genomic.gff.sorted.noregion.gff.gz',
          locationType: 'UriLocation',
        },
        index: {
          location: {
            uri: '/data/GCF_000006765.1/genomic.gff.sorted.noregion.gff.gz.tbi',
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
          refName: 'NC_002516.2',
          start: 0,
          end: 6264404,
          assemblyName: 'GCF_000006765.1',
        },
      ],
      init: {
        assembly: 'GCF_000006765.1',
        loc: 'NC_002516.2:1..5,000',
        tracks: ['RefSeqTrack', 'GFF3GeneTrack'],
      },
    },
  },

  aggregateTextSearchAdapters: [
    {
      type: 'TrixTextSearchAdapter',
      textSearchAdapterId: 'lesb58-index',
      ixFilePath: {
        uri: '/data/GCF_000006765.1/trix/GCF_000006765.1.ix',
        locationType: 'UriLocation',
      },
      ixxFilePath: {
        uri: '/data/GCF_000006765.1/trix/GCF_000006765.1.ixx',
        locationType: 'UriLocation',
      },
      metaFilePath: {
        uri: '/data/GCF_000006765.1/trix/GCF_000006765.1_meta.json',
        locationType: 'UriLocation',
      },
      assemblyNames: ['GCF_000006765.1'],
    },
  ],
};

export default config;
