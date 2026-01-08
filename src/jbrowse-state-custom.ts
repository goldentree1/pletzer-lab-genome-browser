import type { JBrowseConfig, ViewModel } from './types';
import { applyPatch } from 'mobx-state-tree';
import { createViewState } from '@jbrowse/react-linear-genome-view2';

export class GenomeBrowserState {
  private static jbrowseCustomisations: Omit<
    JBrowseConfig,
    'assembly' | 'tracks'
  > = {
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

  private viewState: ViewModel;
  constructor(config: JBrowseConfig) {
    // overwrite conf with customisations
    const newConf = {
      ...config,
      ...GenomeBrowserState.jbrowseCustomisations,
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
    this.viewState = createViewState({ ...newConf, plugins: [] });

    const view = this.viewState.session.views[0];
    if (!view || !view.tracks) return;
    for (const viewTrack of view.tracks) {
      // multiwiggle tracks should have overlapping lines
      for (const display of viewTrack.displays || []) {
        if (display.type === 'MultiLinearWiggleDisplay') {
          applyPatch(display, [
            {
              op: 'replace',
              path: '/rendererTypeNameState',
              value: 'multiline',
            },
          ]);
        }
      }
    }
  }

  public getState() {
    return this.viewState;
  }
}
