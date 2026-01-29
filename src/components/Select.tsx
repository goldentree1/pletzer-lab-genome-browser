interface SelectProps {
  label: string;
  value: string;
  setValue: (value: string) => void;
  options: string[];
  id: string;
  className: string;
}

function Select({ label, value, setValue, options, id }: SelectProps) {
  return (
    <div className="select">
      <label htmlFor={id}>{label}</label>
      <select id={id} value={value} onChange={e => setValue(e.target.value)}>
        {options.map((option, index) => (
          <option key={index} value={option}>
            {option}
          </option>
        ))}
      </select>
    </div>
  );
}

export default Select;
