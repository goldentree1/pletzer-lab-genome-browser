import type { ConfigBuilderOpts, JBrowseConfig, ViewModel } from './types';
import type PluginManager from '@jbrowse/core/PluginManager';
import { applyPatch } from 'mobx-state-tree';
import { createViewState } from '@jbrowse/react-linear-genome-view2';
import Plugin from '@jbrowse/core/Plugin';

/** My custom plugin for JBrowse (allows more customisation) */
export default class MyJbrowsePlugin extends Plugin {
  name = 'MyJbrowsePlugin';

  install(pluginManager: PluginManager) {
    console.log('[MyJbrowsePlugin] installed');

    // example: remove track menu for GFF3 FeatureTracks
    pluginManager.addToExtensionPoint('TrackMenuItems', (items, ctx) => {
      const model = ctx?.model;

      if (
        /** @ts-expect-error fk typescript */
        model?.configuration?.trackType === 'FeatureTrack' &&
        /** @ts-expect-error fk typescript */
        model?.adapterConfig?.type === 'Gff3TabixAdapter'
      ) {
        return []; // meant to remove ... menu but doesnt work
      }

      pluginManager.addToExtensionPoint('FeatureMenuItems', (items, ctx) => {
        /** @ts-expect-error fk typescript */
        if (ctx?.track?.adapterConfig?.type === 'Gff3TabixAdapter') {
          return [];
        }
        return items;
      });

      return items;
    });
  }
}

/** JBrowse config with my custom styling + plugins */
const jbrowseCustomisations: Omit<JBrowseConfig, 'assembly' | 'tracks'> = {
  configuration: {
    // custom styling
    theme: {
      palette: {
        primary: {
          main: '#155B8E',
        },
        secondary: {
          main: '#155B8E',
        },
        tertiary: {
          main: '#e58703',
        },
        quaternary: {
          main: '#22cb04',
        },
      },
      typography: { fontSize: 12 },
      spacing: 5,
      logoPath: {
        uri: '',
      },
    },

    // allow multi-threading (worker)
    rpc: {
      defaultDriver: 'WebWorkerRpcDriver',
    },
  },

  // allow multi-threading (worker)
  makeWorkerInstance: () => {
    return new Worker(new URL('./rpcWorker', import.meta.url), {
      type: 'module',
    });
  },

  // add plugins
  plugins: [MyJbrowsePlugin],
};

/** JBrowse's createViewState() ... but with my customisations */
export function myCreateViewState(config: JBrowseConfig): ViewModel {
  // overwrite conf with customisations
  const newConf = {
    ...config,
    ...jbrowseCustomisations,
  };

  // track patches
  for (const track of newConf.tracks || []) {
    // add reverse/forward strand colouring to genes (GFF3 track)
    if (
      track.type === 'FeatureTrack' &&
      track?.adapter?.type === 'Gff3TabixAdapter'
    ) {
      track.displays = [
        {
          type: 'LinearBasicDisplay',
          displayId: 'GFF3GeneTrack-LinearBasicDisplay',
          height: 130,
          renderer: {
            type: 'SvgFeatureRenderer',
            color1: "jexl:get(feature,'strand')>0?'#00FF00':'#ff0000'",
            color2: "jexl:get(feature,'strand')>0?'#029f02':'#af0101'",
            height: 15,
          },
        },
      ];
    }
  }

  // edit default session settings
  if (newConf.defaultSession?.view) {
    newConf.defaultSession.view.trackLabels = 'overlapping';
    newConf.defaultSession.view.colorByCDS = false;
  }

  // build the views
  const state = createViewState({ ...newConf });

  const view = state.session.views[0];
  if (!view || !view.tracks) return state;
  for (const viewTrack of view.tracks) {
    // multiwiggle tracks should have overlapping lines
    for (const display of viewTrack.displays || []) {
      if (display.type === 'MultiLinearWiggleDisplay') {
        applyPatch(display, [
          { op: 'replace', path: '/rendererTypeNameState', value: 'multiline' },
        ]);
      }
    }
  }

  return state satisfies ViewModel;
}

/** Build a JBrowse-compatible config from my simple config builder */
export function buildConfig(
  {
    firstRegion,
    ncbiName,
    dataDir,
    data: { refSeq, genomic, coverage },
    extras,
  }: ConfigBuilderOpts,
  {
    loc = [0, 5000],
    logScale = true,
    conditionA = [0, 0],
    conditionB = [1, 0],
  } = {},
): JBrowseConfig {
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
          loc: `${firstRegion}:${loc[0] + 1}..${loc[1]}`,
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
  } satisfies JBrowseConfig;

  // add multiwig coverage if at least 2 conditions exist
  if (coverage.length >= 2) {
    conf.tracks.push({
      type: 'MultiQuantitativeTrack',
      trackId: 'multiwig-coverage',
      name: 'Coverage',
      assemblyNames: ['asm'],
      category: ['Coverage'],
      adapter: {
        type: 'MultiWiggleAdapter',
        /** @ts-expect-error works at runtime */
        bigWigs: [
          coverage[conditionA[0]][conditionA[1]],
          coverage[conditionB[0]][conditionB[1]],
        ].map(fname => `${dataDir}/${fname}`),
      },
      displays: [
        {
          type: 'MultiLinearWiggleDisplay',
          displayId: 'Coverage_multiwiggle-MultiLinearWiggleDisplay',
          renderer: { type: 'XYPlotRenderer' },
          height: 300,
          scaleType: logScale ? 'log' : 'linear',
          autoscale: 'global',
        },
      ],
    });

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
