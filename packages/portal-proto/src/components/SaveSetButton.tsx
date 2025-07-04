import React, { useEffect, useState } from "react";
import { TypedUseMutation } from "@reduxjs/toolkit/query/react";
import {
  selectSetsByType,
  useCoreDispatch,
  useCoreSelector,
  addSet,
  SetTypes,
  hideModal,
  CreateSetValueArgs,
  showModal,
  Modals,
} from "@gff/core";
import { showNotification } from "@mantine/notifications";
import { SaveOrCreateEntityModal } from "@/components/Modals/SaveOrCreateEntityModal";
import DarkFunctionButton from "@/components/StyledComponents/DarkFunctionButton";
import { Loader } from "@mantine/core";
interface SaveSetButttonProps {
  readonly disabled: boolean;
  readonly ids: string[];
  readonly hooks: {
    createSet: TypedUseMutation<any, CreateSetValueArgs, any>;
  };
  readonly setType: SetTypes;
  readonly buttonText?: string;
  readonly dismissModal?: boolean;
}

const SaveSetButton: React.FC<SaveSetButttonProps> = ({
  disabled,
  ids,
  hooks,
  setType,
  buttonText = "Save Set",
  dismissModal = false,
}: SaveSetButttonProps) => {
  const dispatch = useCoreDispatch();
  const [showSaveModal, setShowSaveModal] = useState(false);
  const [createSet, response] = hooks.createSet();
  const [setName, setSetName] = useState(null);
  const sets = useCoreSelector((state) => selectSetsByType(state, setType));

  useEffect(() => {
    if (response.isSuccess && setName) {
      dispatch(addSet({ setType, setName, setId: response.data }));
      showNotification({
        message: "Set has been saved.",
        closeButtonProps: { "aria-label": "Close notification" },
      });
      if (dismissModal) {
        dispatch(hideModal());
      }
      setSetName(null);
    } else if (response.isError) {
      dispatch(showModal({ modal: Modals.SaveSetErrorModal }));
    }
  }, [
    response.isSuccess,
    response.isError,
    response.data,
    setName,
    dispatch,
    setType,
    dismissModal,
  ]);

  return (
    <>
      <SaveOrCreateEntityModal
        entity="set"
        initialName=""
        opened={showSaveModal}
        onClose={() => setShowSaveModal(false)}
        onActionClick={(name: string) => {
          setSetName(name);
          createSet({ values: ids, intent: "user", set_type: "mutable" });
        }}
        onNameChange={(name) => !Object.values(sets).includes(name)}
        additionalDuplicateMessage="This will overwrite it."
      />
      <DarkFunctionButton
        data-testid="button-save-set"
        disabled={disabled}
        onClick={() => setShowSaveModal(true)}
        leftSection={
          response?.isLoading ? <Loader size="sm" color="white" /> : undefined
        }
      >
        {buttonText}
      </DarkFunctionButton>
    </>
  );
};

export const SubmitSaveSetButton: React.FC<SaveSetButttonProps> = (
  props: SaveSetButttonProps,
) => (
  <SaveSetButton
    {...props}
    data-testid="button-submit"
    buttonText="Submit"
    dismissModal
  />
);

export default SaveSetButton;
