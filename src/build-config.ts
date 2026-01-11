import type { JBrowseConfig, ConfigBuilderOpts } from './types';

export function buildConfig({
  firstRegion,
  ncbiName,
  dataDir,
  data: { refSeq, genomic, coverage },
  extras,
}: ConfigBuilderOpts): JBrowseConfig {
  if (!dataDir) dataDir = `/data/${ncbiName}`;

  // config to be returned
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
          tracks: ['refseq', 'genomic'],
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

  // monkey patches

  if (coverage.length >= 2) {
    // add multiwig coverage if exists
    conf.tracks.push({
      type: 'MultiQuantitativeTrack',
      trackId: 'multiwig-coverage',
      name: 'Coverage',
      assemblyNames: ['asm'],
      category: ['Coverage'],
      adapter: {
        type: 'MultiWiggleAdapter',
        /** @ts-expect-error idk why type issue, works fine */
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

    // add coverage track to default session
    conf.defaultSession.view.init.tracks.push('multiwig-coverage');
  }

  if (extras) {
    for (const track of extras) {
      conf.tracks.push(track);
      conf.defaultSession.view.init.tracks.push(track.trackId);
    }
  }

  return conf;
}
