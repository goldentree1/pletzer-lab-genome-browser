import type { JBrowseConfig, ViewModel } from './types';
import { applyPatch } from 'mobx-state-tree';
import { createViewState } from '@jbrowse/react-linear-genome-view2';

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
  plugins: [],
};

/**
 * JBrowse's createViewState() ... but with my customisations.
 */
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
          renderer: {
            type: 'SvgFeatureRenderer',
            color1: "jexl:get(feature,'strand')>0?'#00FF00':'#ff0000'",
            color2: "jexl:get(feature,'strand')>0?'#029f02':'#af0101'",
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
  const state = createViewState({ ...newConf, plugins: [] });

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
