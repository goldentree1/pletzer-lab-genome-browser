import type { JBrowseConfig } from './types';

interface Gen {
  ncbiName: string; // e.g., 'GCF_000014625'
  dataDir?: string; // default '/data/$ncbiName/'
  firstRegion: string; // e.g., 'NC_008463.1'
  data: {
    refSeq: string; // default 'REFSEQ.faa.gz'
    genomic: string; // default 'GENOME.gff'
    coverage: string[];
  };
}

export function gen({
  firstRegion,
  ncbiName,
  dataDir,
  data: { refSeq, genomic, coverage },
}: Gen): JBrowseConfig {
  if (!dataDir) dataDir = `/data/${ncbiName}`;

  const conf = {
    assembly: {
      name: 'asm',
      sequence: {
        type: 'ReferenceSequenceTrack',
        trackId: 'refseq',
        adapter: {
          type: 'BgzipFastaAdapter',
          uri: `${dataDir}/${refSeq}`,
        },
      },
    },
    tracks: [
      {
        type: 'FeatureTrack',
        trackId: 'genomic',
        name: 'GFF3 Track',
        assemblyNames: ['asm'],
        adapter: {
          type: 'Gff3TabixAdapter',
          gffGzLocation: {
            uri: `${dataDir}/${genomic}`,
            locationType: 'UriLocation',
          },
          index: {
            location: {
              uri: `${dataDir}/${genomic}.tbi`,
              locationType: 'UriLocation',
            },
          },
        },
      },
    ],
    defaultSession: {
      name: 'default-session',
      margin: 0,
      view: {
        id: 'lgv',
        type: 'LinearGenomeView',
        displayedRegions: [],
        init: {
          assembly: 'asm',
          loc: `${firstRegion}:1..5,000`,
          tracks: ['refseq', 'genomic', 'coverage'],
        },
      },
    },

    aggregateTextSearchAdapters: [
      {
        type: 'TrixTextSearchAdapter',
        textSearchAdapterId: 'text-search',
        ixFilePath: {
          uri: `${dataDir}/trix/${refSeq}.ix`,
          locationType: 'UriLocation',
        },
        ixxFilePath: {
          uri: `${dataDir}/trix/${refSeq}.ixx`,
          locationType: 'UriLocation',
        },
        metaFilePath: {
          uri: `${dataDir}/trix/${refSeq}_meta.json`,
          locationType: 'UriLocation',
        },
        assemblyNames: ['asm'],
      },
    ],
  };

  if (coverage.length > 0) {
    conf.tracks.push({
      type: 'MultiQuantitativeTrack',
      trackId: 'coverage',
      name: 'Coverage',
      assemblyNames: ['asm'],
      category: ['Coverage'],
      adapter: {
        type: 'MultiWiggleAdapter',
        /** @ts-expect-error idk why but moving it out gives bad type */
        bigWigs: coverage.map(c => `${dataDir}/${c}`),
      },
      displays: [
        {
          type: 'MultiLinearWiggleDisplay',
          displayId: 'Coverage_multiwiggle-MultiLinearWiggleDisplay',
          renderer: {
            type: 'XYPlotRenderer',
          },
          scaleType: 'log',
          autoscale: 'global',
        },
      ],
    });
  }

  return conf;
}
