import type { createViewState } from '@jbrowse/react-linear-genome-view2';
type JBrowseConfig = Parameters<typeof createViewState>[0];

const config: JBrowseConfig = {
  assembly: {
    name: 'LESB58',
    sequence: {
      type: 'ReferenceSequenceTrack',
      trackId: 'RefSeqTrack',
      adapter: {
        type: 'BgzipFastaAdapter',
        uri: '/data/NC011770/LESB58_ASM2664v1.fasta.gz',
      },
    },
  },
  tracks: [
    {
      type: 'FeatureTrack',
      trackId: 'GFF3GeneTrack',
      name: 'GFF3 Track',
      assemblyNames: ['LESB58'],
      adapter: {
        type: 'Gff3TabixAdapter',
        gffGzLocation: {
          uri: '/data/NC011770/LESB58_ASM2664v1.gff.sorted.gff.noregion.gz',
          locationType: 'UriLocation',
        },
        index: {
          location: {
            uri: '/data/NC011770/LESB58_ASM2664v1.gff.sorted.gff.noregion.gz.tbi',
            locationType: 'UriLocation',
          },
        },
      },
      displays: [
        {
          type: 'LinearBasicDisplay',
          displayId: 'GFF3GeneTrack-LinearBasicDisplay',
          renderer: {
            type: 'SvgFeatureRenderer',
            color1: "jexl:get(feature,'strand')>0?'#00FF00':'red'",
          },
        },
      ],
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
          refName: 'NC_011770.1',
          start: 0,
          end: 6601757,
          assemblyName: 'LESB58',
        },
      ],
      init: {
        assembly: 'LESB58',
        loc: 'NC_011770.1:1..5,000',
        tracks: ['RefSeqTrack', 'GFF3GeneTrack'],
      },
    },
  },

  aggregateTextSearchAdapters: [
    {
      type: 'TrixTextSearchAdapter',
      textSearchAdapterId: 'lesb58-index',
      ixFilePath: {
        uri: '/data/NC011770/trix/assembly1.ix',
        locationType: 'UriLocation',
      },
      ixxFilePath: {
        uri: '/data/NC011770/trix/assembly1.ixx',
        locationType: 'UriLocation',
      },
      metaFilePath: {
        uri: '/data/NC011770/trix/assembly1_meta.json',
        locationType: 'UriLocation',
      },
      assemblyNames: ['LESB58'],
    },
  ],

  configuration: {
    rpc: {
      defaultDriver: 'WebWorkerRpcDriver',
    },
  },

  makeWorkerInstance: () => {
    return new Worker(new URL('./rpcWorker', import.meta.url), {
      type: 'module',
    });
  },
};

export default config;
