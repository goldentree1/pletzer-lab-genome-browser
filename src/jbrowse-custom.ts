import type {
  JBrowseConfig,
  ViewModel,
  JBrowseCustomConfigHybrid,
} from './types';
// import type PluginManager from '@jbrowse/core/PluginManager';
import { applyPatch } from 'mobx-state-tree';
import { createViewState } from '@jbrowse/react-linear-genome-view2';

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
    // plugins: [MyJbrowsePlugin],
  };

/** JBrowse's createViewState() ... but with my customisations */
export function myCreateViewState(config: JBrowseConfig): ViewModel {
  const newConf: JBrowseConfig = {
    ...config,
    ...staticJBrowseCustomisations,
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
        // // 2. Use the built-in JBrowse action to change height
        // if (typeof display.setHeight === 'function') {
        //   display.setHeight(200); // Set your desired pixel height here
        // } else {
        //   // Fallback for some display types
        //   display.height = 200;
        // }
      }
    }
  }

  return state satisfies ViewModel;
}

/** Build a JBrowse-compatible config from my simple config builder */
export function buildConfig(
  {
    firstRegion,
    genomeName,
    dataDir,
    trixName,
    norms,
    genesLabelTypes,
    data: { refSeq, genomic, experiments },
    extras,
  }: JBrowseCustomConfigHybrid,
  {
    loc = [0, 5000],
    logScale = true,
    experiment = Object.keys(experiments)[0],
    conditionA = [0, 0],
    conditionB = [1, 0],
    extraConditions = [],
    normType = 'none',
    genesLabelType = 'name',
  },
): JBrowseConfig {
  if (!trixName) trixName = genomeName;
  const baseUri = new URL('.', window.location.href).href;
  if (!dataDir) dataDir = `${baseUri}data/${genomeName}`;
  if (!refSeq) refSeq = 'refseq.fna.gz';

  const conf = {
    assembly: {
      name: 'asm',
      sequence: {
        type: 'ReferenceSequenceTrack',
        trackId: 'refseq',
        name: `Reference Sequence`,
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
          type: 'Gff3Adapter',
          gffLocation: {
            uri: `${dataDir}/${genomic}`, // plain .gff
            locationType: 'UriLocation',
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

  const exp = experiments[experiment];
  if (!exp) {
    throw new Error(`Unknown experiment: ${experiment}`);
  }

  const {
    coverage,
    // coverage_condition_names
  } = exp;

  // add multiwig coverage if at least 2 conditions exist
  if (coverage.length) {
    const trackId = `multiwig-coverage-${experiment}`;
    const all = [
      coverage[conditionA[0]]?.[conditionA[1]],
      coverage[conditionB[0]]?.[conditionB[1]],
      ...((extraConditions as [number, number][]) ?? []).map(
        ([i, j]) => coverage[i]?.[j],
      ),
    ].filter(Boolean);

    conf.tracks.push({
      type: 'MultiQuantitativeTrack',
      trackId,
      name: `Coverage (${experiment})`,
      assemblyNames: ['asm'],
      category: ['Coverage', experiment],
      adapter: {
        type: 'MultiWiggleAdapter',
        /** @ts-expect-error works at runtime */
        bigWigs: all.map(fname => {
          const fpath = `${dataDir}/${fname}`;
          if (normType.toLowerCase() === 'cpm' && norms.includes('cpm')) {
            return `${fpath.substring(0, fpath.length - 3)}.cpm.bw`;
          }
          return fpath;
        }),
      },

      displays: [
        {
          type: 'MultiLinearWiggleDisplay',
          displayId: `${trackId}-MultiLinearWiggleDisplay`,
          renderer: { type: 'XYPlotRenderer', height: 600 },
          height: 380,
          scaleType: logScale ? 'log' : 'linear',
          autoscale: 'global',
        },
      ],
    });

    conf.defaultSession.view.init.tracks.push(trackId);
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
    const finalLbl = genesLabelTypes?.includes(genesLabelType)
      ? genesLabelType
      : 'name';
    if (
      track.type === 'FeatureTrack' &&
      track?.adapter?.type === 'Gff3Adapter'
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
