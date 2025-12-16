import type { createViewState } from '@jbrowse/react-linear-genome-view2';
type JBrowseConfig = Parameters<typeof createViewState>[0];

const config: JBrowseConfig = {
  assembly: {
    name: 'NC008463',
    sequence: {
      type: 'ReferenceSequenceTrack',
      trackId: 'RefSeqTrack',
      adapter: {
        type: 'BgzipFastaAdapter',
        uri: '/data/NC008463/GCF_000014625.1_ASM1462v1_genomic.fna.gz',
      },
    },
  },
  tracks: [
    {
      type: 'FeatureTrack',
      trackId: 'GFF3GeneTrack',
      name: 'GFF3 Track',
      assemblyNames: ['NC008463'],
      adapter: {
        type: 'Gff3TabixAdapter',
        gffGzLocation: {
          uri: '/data/NC008463/genomic.gff.sorted.noregion.gz',
          locationType: 'UriLocation',
        },
        index: {
          location: {
            uri: '/data/NC008463/genomic.gff.sorted.noregion.gz.tbi',
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
          refName: 'NC_008463.1',
          start: 0,
          end: 6537648,
          assemblyName: 'NC008463',
        },
      ],
      init: {
        assembly: 'NC008463',
        loc: 'NC_008463.1:1..5,000',
        tracks: ['RefSeqTrack', 'GFF3GeneTrack'],
      },
    },
  },

  // aggregateTextSearchAdapters: [
  //   {
  //     type: 'TrixTextSearchAdapter',
  //     textSearchAdapterId: 'NC008463-index',
  //     ixFilePath: {
  //       uri: '/data/NC011770/trix/assembly1.ix',
  //       locationType: 'UriLocation',
  //     },
  //     ixxFilePath: {
  //       uri: '/data/NC011770/trix/assembly1.ixx',
  //       locationType: 'UriLocation',
  //     },
  //     metaFilePath: {
  //       uri: '/data/NC011770/trix/assembly1_meta.json',
  //       locationType: 'UriLocation',
  //     },
  //     assemblyNames: ['NC008463'],
  //   },
  // ],

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
