import type { ViewModel } from './types';
import { useState, useEffect, useRef } from 'react';
import { JBrowseLinearGenomeView } from '@jbrowse/react-linear-genome-view2';
// @ts-expect-error no types for font
import '@fontsource/roboto';
import myConf from './config';
import { buildConfig, myCreateViewState } from './jbrowse-custom';
import { useStoredStateBoolean, useStoredStateString } from './utils';
// import Select from './components/Select';
import { autorun } from 'mobx';

// eslint-disable-next-line @typescript-eslint/no-explicit-any
function getLocString(linearView: any) {
  if (!linearView.displayedRegions.length) return;
  const region = linearView.displayedRegions[0];

  const bpPerPx = linearView.bpPerPx;
  const offsetPx = linearView.offsetPx;
  const widthPx = linearView.width;

  const start = region.start + Math.floor(offsetPx * bpPerPx) + 1;
  const end = region.start + Math.floor((offsetPx + widthPx) * bpPerPx);

  const s = Math.max(start, 1);
  const e = Math.min(end, region.end);

  return region.reversed
    ? `${region.refName}:${e}..${s}`
    : `${region.refName}:${s}..${e}`;
}

function App() {
  const genomes: string[] = Object.keys(myConf).sort();
  const [genome, setGenome] = useStoredStateString(
    'pletzer-lab-genome-browser:genome',
    genomes[0],
  );

  const [viewState, setViewState] = useState<ViewModel>();
  const [conditionA, setConditionA] = useState<[number, number]>([0, 0]);
  const [conditionB, setConditionB] = useState<[number, number]>([1, 0]);
  const [norms, setNorms] = useState(['none']);
  const loc = useRef<string | null>(null);

  const [genesLabelType, setGenesLabelType] = useStoredStateString(
    'pletzer-lab-genome-browser:genesLabelType',
    'name',
  );
  const [normType, setNormType] = useStoredStateString(
    'pletzer-lab-genome-browser:normType',
    'none',
  );
  const [logScaling, setLogScaling] = useStoredStateBoolean(
    'pletzer-lab-genome-browser:logScaling',
    true,
  );
  const [globalScaling, setGlobalScaling] = useStoredStateBoolean(
    'pletzer-lab-genome-browser:globalScaling',
    true,
  );
  const [colorByCds, setColorByCDS] = useStoredStateBoolean(
    'pletzer-lab-genome-browser:colorByCDS',
    false,
  );

  // DOM manipulation works for some JBrowse features - like scaling and CDS colouring.
  // Smoother and faster than a full-refresh, so use it where possible.
  useEffect(() => {
    const linearView = viewState?.session.views.find(
      v => v.type === 'LinearGenomeView',
    );
    if (!linearView) return;

    linearView.tracks.forEach(track => {
      /** @ts-expect-error display is 'any' */
      track.displays.forEach(display => {
        if (display.type.includes('Wiggle')) {
          const wantedScale = logScaling ? 'log' : 'linear';
          if (display.scaleType !== wantedScale)
            display.setScaleType(wantedScale);
          if (display.autoscale !== (globalScaling ? 'global' : 'local'))
            display.setAutoscale(globalScaling ? 'global' : 'local');
        }
      });
    });

    linearView.setColorByCDS(colorByCds);
  }, [viewState, logScaling, globalScaling, colorByCds]);

  useEffect(() => {
    const linearView = viewState?.session.views.find(
      v => v.type === 'LinearGenomeView',
    );
    if (!linearView) return;

    const dispose = autorun(() => {
      loc.current = getLocString(linearView) || null;
    });

    return () => dispose();
  }, [viewState]);

  // revert default conditions + loc on genome change
  useEffect(() => {
    setConditionA([0, 0]);
    setConditionB([1, 0]);
    loc.current = null;
  }, [genome]);

  // full-refresh jbrowse required for genome or conditions changes
  useEffect(() => {
    const newLoc = [0, 5000];
    if (loc.current) {
      const [newStartLoc, newEndLoc] = loc.current
        .split(':')[1]
        .split('..')
        .map(Number);
      newLoc[0] = newStartLoc;
      newLoc[1] = newEndLoc;
      loc.current = null;
    }

    const config = myConf[genome];
    if (!config) {
      setGenome(genomes[0]);
      return;
    }

    setNorms(config.norms);

    if (!config.norms.includes(normType)) {
      setNormType(config.norms[0]);
      return;
    }

    console.log('using loc:', newLoc);

    const state = myCreateViewState(
      buildConfig(config, {
        conditionA,
        conditionB,
        loc: newLoc, // TODO: keep loc constant for condition changes only
        normType,
        genesLabelType,
      }),
    );

    setViewState(state);
  }, [genome, conditionA, conditionB, normType, genesLabelType]); // eslint-disable-line react-hooks/exhaustive-deps

  const coverage = viewState ? myConf[genome].data.coverage : [];
  const coverageConditionNames = viewState
    ? myConf[genome].data.coverage_condition_names
    : [];

  return (
    <>
      <header className="header">
        <nav className="header-genome-chooser">
          <div className="header-title">
            <img src="./pletzerlab-icon.webp" alt="Pletzer Lab Icon" />
            {/*<h1>Pletzer Lab Genome Browser</h1>*/}
          </div>
          <GenomeSelector
            genomes={genomes}
            selected={genome}
            onChange={setGenome}
          />
          {coverage.length >= 1 && (
            <div className="header-condition-chooser">
              <ConditionsSelect
                id="select-condition-a"
                label="Condition"
                coverage={coverage}
                coverageConditionNames={coverageConditionNames}
                value={conditionA}
                onChange={setConditionA}
              />
              {coverage.length > 1 && <span>vs.</span>}
              {coverage.length > 1 && (
                <ConditionsSelect
                  id="select-condition-b"
                  label="Condition"
                  coverage={coverage}
                  coverageConditionNames={coverageConditionNames}
                  value={conditionB}
                  onChange={setConditionB}
                />
              )}
            </div>
          )}
        </nav>

        <div className="header-buttons-container">
          <div className="header-buttons">
            <GenesLabelTypeSelect
              geneLabelTypes={['name', 'locus_tag', 'old_locus_tag']}
              geneLabelType={genesLabelType}
              setGeneLabelType={setGenesLabelType}
            />
            <NormTypeSelect
              norms={norms}
              normType={normType}
              setNormType={setNormType}
            />
            <BooleanCheckbox
              label="Log2"
              checked={logScaling}
              onChange={setLogScaling}
              hoverDescription="Use base-2 log scale"
            />
            <BooleanCheckbox
              label="Global Y-axis"
              checked={globalScaling}
              onChange={setGlobalScaling}
              hoverDescription="Fix Y-axis between global minimum and maximum values"
            />
            <BooleanCheckbox
              label="CDS+Amino"
              checked={colorByCds}
              onChange={setColorByCDS}
              hoverDescription="Colour by CDS and amino acids (shows codon alignments by colour)"
            />
          </div>
        </div>
      </header>

      <div className="jbrowse-container">
        {viewState && <JBrowseLinearGenomeView viewState={viewState} />}
      </div>
    </>
  );
}

export default App;

function GenomeSelector({
  genomes,
  selected,
  onChange,
}: {
  genomes: string[];
  selected: string;
  onChange: (b: string) => void;
}) {
  return (
    <select value={selected} onChange={e => onChange(e.target.value)}>
      {genomes.map(b => (
        <option key={b} value={b}>
          {b}
        </option>
      ))}
    </select>
  );
}

function ConditionsSelect({
  coverage,
  coverageConditionNames,
  value,
  onChange,
  id,
}: {
  label: string;
  coverage: string[][];
  coverageConditionNames: string[];
  value: [number, number];
  onChange: (val: [number, number]) => void;
  className?: string;
  id: string;
}) {
  return (
    <select
      className="condition-select"
      id={id}
      value={value.join(',')}
      onChange={e => {
        const [c, r] = e.target.value.split(',').map(Number);
        onChange([c, r]);
      }}
    >
      {coverage.map((arr, i) => (
        <optgroup
          key={`${id}--group-${i}`}
          label={`Condition ${i + 1}: ${coverageConditionNames[i]}`}
        >
          {arr.map((filepath, j) => (
            <option key={`${id}--group-${i}-item-${j}`} value={`${i},${j}`}>
              {optname(filepath)}
            </option>
          ))}
        </optgroup>
      ))}
    </select>
  );

  function optname(filepath: string) {
    let filename = filepath.substring(filepath.lastIndexOf('/') + 1);
    if (filename.endsWith('.average.bw')) {
      filename = `${filename.substring(0, filename.length - 11)} (average)`;
    } else if (filename.match(/\d+\.bw$/)) {
      const noBw = filename.substring(0, filename.length - 3);
      const replicateN = parseInt(noBw.substring(noBw.lastIndexOf('.') + 1));
      filename = `${noBw.substring(0, noBw.lastIndexOf('.'))} (replicate #${replicateN})`;
    }
    return `${filename}`;
  }
}

function NormTypeSelect({
  norms,
  normType,
  setNormType,
}: {
  norms: string[];
  normType: string;
  setNormType: (value: string) => void;
}) {
  return (
    <div>
      <label htmlFor="norm-type-select">Normalisation: </label>
      <select
        id="norm-type-select"
        value={normType}
        className="norm-type-select"
        onChange={e => setNormType(e.target.value)}
      >
        {norms.map((norm, index) => (
          <option key={index} value={norm}>
            {norm}
          </option>
        ))}
      </select>
    </div>
  );
}

function GenesLabelTypeSelect({
  geneLabelTypes,
  geneLabelType,
  setGeneLabelType,
}: {
  geneLabelTypes: string[];
  geneLabelType: string;
  setGeneLabelType: (value: string) => void;
}) {
  return (
    <div>
      <label htmlFor="genes-label-type-select">Gene Labels:</label>
      <select
        id="genes-label-type-select"
        value={geneLabelType}
        className="genes-label-type-select"
        onChange={e => setGeneLabelType(e.target.value)}
      >
        {geneLabelTypes.map((geneLabelType, index) => (
          <option key={index} value={geneLabelType}>
            {geneLabelType}
          </option>
        ))}
      </select>
    </div>
  );
}

function BooleanCheckbox({
  label,
  checked,
  onChange,
  hoverDescription,
}: {
  label: string;
  checked: boolean;
  onChange: (b: boolean) => void;
  hoverDescription: string;
}) {
  const id = label.replace(/\s+/g, '').toLowerCase();
  return (
    <div className="checkbox">
      <label htmlFor={id} title={hoverDescription}>
        {label}
      </label>
      <input
        type="checkbox"
        id={id}
        checked={checked}
        onChange={e => onChange(e.target.checked)}
        title={hoverDescription}
      />
    </div>
  );
}
