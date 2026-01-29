import type { JBrowseCustomConfig, JBrowseConfig, ViewModel } from './types';
import type PluginManager from '@jbrowse/core/PluginManager';
import { applyPatch } from 'mobx-state-tree';
import { createViewState } from '@jbrowse/react-linear-genome-view2';
import Plugin from '@jbrowse/core/Plugin';
/**@ts-ignore */
import TestPlugin from './plugin-test';

// /** My custom plugin for JBrowse (allows more customisation) */
/**@ts-ignore */
import { getSnapshot } from 'mobx-state-tree';
export default class MyJbrowsePlugin extends Plugin {
  name = 'MyJbrowsePlugin';
  /**@ts-ignore */
  install(pluginManager: PluginManager) {
    // Keep this just to see it fire
    console.log('--- PLUGIN INSTALLING ---');
  }

  configure(pluginManager: PluginManager) {
    console.log('--- PLUGIN CONFIGURING (System should be ready) ---');

    const epNames = Object.keys(pluginManager);

    console.log('EPs available now:', epNames);
    // console.log(JSON.stringify(pluginManager));

    // Now try to add the menu item
    // Usually 'TrackMenuItems' or 'LinearGenomeView-TrackMenu'
    // const menuName = epNames.includes('LinearGenomeView-TrackMenu')
    //   ? 'LinearGenomeView-TrackMenu'
    //   : 'TrackMenuItems';

    // if (epNames.length > 0) {
    //   pluginManager.addToExtensionPoint(menuName, (items, ctx) => {
    //     console.log('Menu hook finally fired!');
    //     return [
    //       /** @ts-expect-error */
    //       ...items,
    //       {
    //         label: 'Force Locus Tags',
    //         onClick: () => {
    //           // Brute force patch
    //           const target = ctx.track || ctx.model;
    //           /** @ts-expect-error */
    //           applyPatch(target.configuration.renderer, [
    //             {
    //               op: 'replace',
    //               path: '/labels/name',
    //               value: "jexl:get(feature, 'old_locus_tag')",
    //             },
    //           ]);
    //         },
    //       },
    //     ];
    //   });
    // }
  }
}
/** JBrowse config with my custom styling + plugins */
const staticJBrowseCustomisations: Omit<JBrowseConfig, 'assembly' | 'tracks'> =
  {
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
  const newConf: JBrowseConfig = {
    ...config,
    ...staticJBrowseCustomisations,
    plugins: [MyJbrowsePlugin],
  };
  // build the views
  const state = createViewState(newConf);

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
    trixName,
    norms,
    genesLabelTypes,
    data: { refSeq, genomic, coverage },
    extras,
  }: JBrowseCustomConfig,
  {
    loc = [0, 5000],
    logScale = true,
    conditionA = [0, 0],
    conditionB = [1, 0],
    normType = 'none',
    genesLabelType = 'name',
  } = {},
): JBrowseConfig {
  if (!trixName) trixName = ncbiName;
  const baseUri = new URL('.', window.location.href).href;
  if (!dataDir) dataDir = `${baseUri}data/${ncbiName}`;
  if (!refSeq) refSeq = 'refseq.fna.gz';

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
        name: 'Genes',
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
          uri: `${dataDir}/trix/${trixName}.ix`,
          locationType: 'UriLocation',
        },
        ixxFilePath: {
          uri: `${dataDir}/trix/${trixName}.ixx`,
          locationType: 'UriLocation',
        },
        metaFilePath: {
          uri: `${dataDir}/trix/${trixName}_meta.json`,
          locationType: 'UriLocation',
        },
        assemblyNames: ['asm'],
      },
    ],
  } satisfies JBrowseConfig; // satisfies so it knows it's not a proper volatile config yet

  // add multiwig coverage if at least 2 conditions exist
  if (coverage.length) {
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

          // quick hacky fix to get only one coverage condition working
          coverage.length >= 2 ? coverage[conditionB[0]][conditionB[1]] : null,
        ]
          .filter(Boolean)
          .map(fname => {
            const fpath = `${dataDir}/${fname}`;
            if (normType.toLowerCase() === 'cpm' && norms.includes('cpm')) {
              return `${fpath?.substring(0, fpath.length - 3)}.cpm.bw`;
            }
            return `${dataDir}/${fname}`;
          }),
      },
      displays: [
        {
          type: 'MultiLinearWiggleDisplay',
          displayId: 'Coverage_multiwiggle-MultiLinearWiggleDisplay',
          renderer: { type: 'XYPlotRenderer' },
          height: 340,
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

  const finalConf = conf as JBrowseConfig; // as, so we can add new properties
  // track patches
  for (const track of finalConf.tracks || []) {
    // add reverse/forward strand colouring to genes (GFF3 track)
    console.log('requested-lbl', genesLabelType);
    const finalLbl = genesLabelTypes?.includes(genesLabelType)
      ? genesLabelType
      : 'name';
    console.log('finallbl', finalLbl);
    if (
      track.type === 'FeatureTrack' &&
      track?.adapter?.type === 'Gff3TabixAdapter'
    ) {
      track.displays = [
        {
          type: 'LinearBasicDisplay',
          displayId: 'GFF3GeneTrack-LinearBasicDisplay',
          height: 175,
          renderer: {
            type: 'SvgFeatureRenderer',
            color1: "jexl:get(feature,'strand')>0?'#00FF00':'#ff0000'",
            color2: "jexl:get(feature,'strand')>0?'#029f02':'#af0101'",
            height: 15,
            labels: {
              // main label
              name: `jexl:get(feature, '${finalLbl}')`,
              // description: "jexl:get(feature, 'old_locus_tag')", // secondary label (we may want?)
            },
          },
        },
      ];
    }
  }

  // edit default session settings
  if (finalConf.defaultSession?.view) {
    finalConf.defaultSession.view.trackLabels = 'overlapping';
    finalConf.defaultSession.view.colorByCDS = false;
  }

  return finalConf;
}
