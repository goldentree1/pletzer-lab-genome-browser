import type { JBrowseConfig, JBrowseConfigBuilderOpts } from './types';

const defaultConf: Omit<JBrowseConfig, 'assembly' | 'tracks'> = {
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
};

/**
 * Transform the given configuration into my nicely customised and styled configuration
 */
export function makeConfig(config: JBrowseConfig): JBrowseConfig {
  // do as much as we can by overwriting with defaults
  const newConf = {
    ...config,
    ...defaultConf,
  };

  // dynamic track patches
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
            color1: "jexl:get(feature,'strand')>0?'#00FF00':'red'",
            color2: "jexl:get(feature,'strand')>0?'#00FF00':'red'",
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

  return newConf;
}

export function makeConfigFromOpts(
  opts: JBrowseConfigBuilderOpts,
): JBrowseConfig {
  console.log(opts);

  return {
    assembly: {
      name: opts.assemblyName,
      sequence: {
        type: 'ReferenceSequenceTrack',
        trackId: 'RefSeqTrack',
        adapter: {
          type: 'BgzipFastaAdapter',
          uri: `/data/${opts.dataDirName}/LESB58_ASM2664v1.fasta.gz`,
        },
      },
    },
    tracks: [],
    ...defaultConf,
  };
}
