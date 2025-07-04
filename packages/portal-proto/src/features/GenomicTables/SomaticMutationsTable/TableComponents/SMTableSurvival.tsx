import { SurvivalChartIcon } from "@/utils/icons";
import ToggledCheck from "../../SharedComponent/ToggledCheck";
import { Survival } from "../types";

const SMTableSurvival = ({
  affectedCasesInCohort,
  survival,
  proteinChange,
  handleSurvivalPlotToggled,
}: {
  affectedCasesInCohort: {
    numerator: number;
    denominator: number;
  };
  survival: Survival;
  proteinChange: {
    symbol: string;
    aaChange: string;
    geneId: string;
  };
  handleSurvivalPlotToggled: (
    symbol: string,
    name: string,
    field: string,
  ) => void;
}): JSX.Element => {
  const { numerator } = affectedCasesInCohort ?? {
    numerator: 0,
  };
  const disabled = numerator < 10;
  const isActive = survival.checked;
  const tooltip = disabled
    ? `Not enough data`
    : isActive
    ? `Remove ${survival.name} from plot`
    : `Plot ${survival.name}`;

  return (
    <ToggledCheck
      ariaText={`Toggle survival plot for ${proteinChange.symbol} ${proteinChange.aaChange} mutation`}
      margin="ml-0.5"
      isActive={survival?.checked}
      icon={<SurvivalChartIcon size={24} aria-hidden="true" />}
      survivalProps={{ plot: "gene.ssm.ssm_id" }}
      selected={survival as unknown as Record<string, string>} // need to fix this
      disabled={disabled}
      handleSwitch={handleSurvivalPlotToggled}
      tooltip={tooltip}
    />
  );
};

export default SMTableSurvival;
