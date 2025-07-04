import { AlertIcon } from "@/utils/icons";

export const WarningBanner = ({ text }: { text: string }): JSX.Element => (
  <div className="flex h-12 border border-warningColor rounded-none">
    <div className="flex h-full w-12 bg-warningColor justify-center items-center">
      <AlertIcon color="white" className="h-6 w-6" aria-label="Warning" />
    </div>
    <div className="bg-[#FFAD0D33] h-full w-full flex items-center pl-4">
      <span data-testid="text-warning-banner" className="text-sm">
        {text}
      </span>
    </div>
  </div>
);
