from pydantic import BaseModel, Field
import json
import pandas as pd


class ProteinRMSD(BaseModel):
    """
    Pairwise protein structural RMSD between two PreppedComplex structures.
    Mirrors the MoleculeSimilarity interface so results can be handled uniformly.
    """

    Reference_Structure: str = Field(
        ..., description="ID of the reference protein structure"
    )
    Query_Structure: str = Field(
        ..., description="ID of the query (mobile) protein structure"
    )
    RMSD: float = Field(
        ..., ge=0, description="Cα RMSD in Ångströms after superposition"
    )
    Binding_Site_Only: bool = Field(
        False,
        description="True if only binding-site residues were used for alignment; "
        "False if the full chain-A was used",
    )

    def save(self, path):
        with open(path, "w") as f:
            json.dump(self.model_dump(), f, indent=4)

    @classmethod
    def load(cls, path):
        with open(path, "r") as f:
            return cls(**json.load(f))

    def __str__(self):
        scope = "binding-site" if self.Binding_Site_Only else "full chain-A"
        return (
            f"Protein RMSD ({scope}) between {self.Reference_Structure} "
            f"and {self.Query_Structure}: {self.RMSD:.4f} Å"
        )

    @classmethod
    def construct_dataframe(cls, rmsd_list: list["ProteinRMSD"]) -> pd.DataFrame:
        return pd.DataFrame.from_records([r.model_dump() for r in rmsd_list])

    @classmethod
    def from_superposition(
        cls,
        ref_id: str,
        mobile_id: str,
        rmsd: float,
        binding_site_only: bool,
    ) -> "ProteinRMSD":
        """Convenience constructor that accepts the raw outputs of superpose_molecule."""
        return cls(
            Reference_Structure=ref_id,
            Query_Structure=mobile_id,
            RMSD=round(rmsd, 4),
            Binding_Site_Only=binding_site_only,
        )
