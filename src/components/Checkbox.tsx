export default function Checkbox({
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
