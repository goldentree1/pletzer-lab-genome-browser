import type { ViewModel } from './types';
import { useState, useEffect } from 'react';
import { JBrowseLinearGenomeView } from '@jbrowse/react-linear-genome-view2';
// @ts-expect-error no types for font
import '@fontsource/roboto';
import myConf from './config';
import { buildConfig, myCreateViewState } from './jbrowse-custom';
import { useStoredStateBoolean, useStoredStateString } from './utils';

function App() {
  const bacteria: string[] = Object.keys(myConf).sort();
  const [bacterium, setBacterium] = useStoredStateString(
    'pletzer-lab-genome-browser:bacterium',
    bacteria[0],
  );

  const [viewState, setViewState] = useState<ViewModel>();
  const [conditionA, setConditionA] = useState<[number, number]>([0, 0]);
  const [conditionB, setConditionB] = useState<[number, number]>([1, 0]);
  const [norms, setNorms] = useState(['none']);

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

  // use default conditions on bacterium change
  useEffect(() => {
    setConditionA([0, 0]);
    setConditionB([1, 0]);
  }, [bacterium]);

  // full-refresh jbrowse required for bacterium or conditions changes
  useEffect(() => {
    const config = myConf[bacterium];
    if (!config) {
      setBacterium(bacteria[0]);
      return;
    }

    setNorms(config.norms);

    if (!config.norms.includes(normType)) {
      setNormType(config.norms[0]);
      return;
    }

    const state = myCreateViewState(
      buildConfig(config, {
        conditionA,
        conditionB,
        loc: [0, 5000], // TODO: keep loc constant for condition changes only
        normType,
      }),
    );

    setViewState(state);
  }, [bacterium, conditionA, conditionB, normType]);

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

  const coverage = viewState ? myConf[bacterium].data.coverage : [];
  const coverageConditionNames = viewState
    ? myConf[bacterium].data.coverage_condition_names
    : [];

  return (
    <>
      <header className="header">
        <div className="header-title">
          <img src="./pletzerlab-icon.webp" alt="Pletzer Lab Icon" />
          <h1>Pletzer Lab Genome Browser</h1>
        </div>

        <nav className="header-genome-chooser">
          <GenomeSelector
            bacteria={bacteria}
            selected={bacterium}
            onChange={setBacterium}
          />

          {coverage.length >= 1 && (
            <div className="header-condition-chooser">
              <ConditionSelect
                id="select-condition-a"
                label="Condition"
                coverage={coverage}
                coverageConditionNames={coverageConditionNames}
                value={conditionA}
                onChange={setConditionA}
              />
              {coverage.length > 1 && <span>vs.</span>}
              {coverage.length > 1 && (
                <ConditionSelect
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
            <div className="norm-type-select">
              <label htmlFor="norm-type">Normalization Type: </label>
              <NormTypeSelect
                norms={norms}
                normType={normType}
                setNormType={setNormType}
              />
            </div>
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
  bacteria,
  selected,
  onChange,
}: {
  bacteria: string[];
  selected: string;
  onChange: (b: string) => void;
}) {
  return (
    <select value={selected} onChange={e => onChange(e.target.value)}>
      {bacteria.map(b => (
        <option key={b} value={b}>
          {b}
        </option>
      ))}
    </select>
  );
}

function ConditionSelect({
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
    } else if (filename.endsWith('.average.cpm.bw')) {
      filename = `${filename.substring(0, filename.length - 15)} (CPM average)`;
    } else if (filename.match(/\d+\.bw$/)) {
      const noBw = filename.substring(0, filename.length - 3);
      const replicateN = parseInt(noBw.substring(noBw.lastIndexOf('.') + 1));
      filename = `${noBw.substring(0, noBw.lastIndexOf('.'))} (replicate #${replicateN})`;
    } else if (filename.match(/\d+\.cpm\.bw$/)) {
      const noBw = filename.substring(0, filename.length - 7);
      const replicateN = parseInt(noBw.substring(noBw.lastIndexOf('.') + 1));
      filename = `${noBw.substring(0, noBw.lastIndexOf('.'))} (CPM replicate #${replicateN})`;
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
    <select value={normType} onChange={e => setNormType(e.target.value)}>
      {norms.map((norm, index) => (
        <option key={index} value={norm}>
          {norm}
        </option>
      ))}
    </select>
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
