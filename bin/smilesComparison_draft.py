import sys
import rdkit
import json
import re
from rdkit import rdBase
from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit import DataStructs
from rdkit.Chem import rdMolDescriptors


# takes string and replaces everything within [] with *, or if no [ exists replaces * with [*]. returns string
def remove_positional_information(x):
    replaced_start = False
    replaced_end = False
    if x[0] == '*':
        x = '[*]' + x[1:]
        replaced_start = True
    if x[(len(x) - 1)] == '*':
        x = x[0:(len(x) - 1)] + '[*]'
        replaced_end = True

    if replaced_start == False and x[0] == '[' and x[1].isdigit():
        for i in [1, 2, 3]:
            if x[(i + 1):(i + 3)] == "*]":
                x = '[*]' + x[(i + 3):]
                replaced_start = True
                break
            elif not x[(1 + i)].isdigit():

                print(x[(i + 1):(i + 3)])
                break

    if replaced_end == False and x[(len(x) - 2):] == '*]' and x[(len(x) - 3)].isdigit():
        for i in [1, 2, 3]:
            if x[(len(x) - 3 - i)] == "[":
                x = x[:(len(x) - 3 - i)] + '[*]'
                replaced_end = True
                break
            elif not x[(len(x) - 3 - i)].isdigit():
                break
    return x


# calculate summary statistics for substituents returned by Chem.MolToSmiles(Chem.ReplaceCore()) and gives position for
# each substituent as dict.
def summarize_substitutians_across_molecule(substituents):
    i = 0
    subs_w_pos = substituents.split(".")

    if (len(subs_w_pos[0]) > 0 and subs_w_pos[0] != ''):
        subs_unified = list(map(remove_positional_information, subs_w_pos))
        sub_dict = {}

        for i in range(len(subs_unified)):
            if subs_unified[i] in sub_dict.keys():
                sub_dict[subs_unified[i]].append(subs_w_pos[i][1:subs_w_pos[i].find("*")])
            else:
                sub_dict[subs_unified[i]] = [subs_w_pos[i][1:subs_w_pos[i].find("*")]]

        return sub_dict
    elif len(subs_w_pos[0]) == 0:
        return ""


def smiles_comparisson(smiles1, smiles2, id1, id2, min_overlapping_atoms=5, MCStimeout=500, fp_maxPath=7,
                       fp_fpSize=2048):
    smiles1 = smiles1.replace('@', '')
    smiles2 = smiles2.replace('@', '')

    # remove implicit and explicit hydrogens and generate formulas
    mol1 = Chem.rdmolops.RemoveHs(Chem.MolFromSmiles(smiles1))
    mol2 = Chem.rdmolops.RemoveHs(Chem.MolFromSmiles(smiles2))

    # calculate overlap between molecules
    overlap = rdFMCS.FindMCS([mol1, mol2],
                             atomCompare=rdFMCS.AtomCompare.CompareElements,
                             bondCompare=rdFMCS.BondCompare.CompareOrderExact,
                             maximizeBonds=True,
                             matchValences=False,
                             timeout=MCStimeout
                             )

    # only accept overlap with more than min_overlapping_atoms atoms
    if (overlap.canceled != True and overlap.numAtoms >= min_overlapping_atoms):

        overlap_dict = {"smartsString": overlap.smartsString,
                        "numAtoms": overlap.numAtoms,
                        "numBonds": overlap.numBonds,
                        "canceled": overlap.canceled}

        # calculate remaining chemical groups after removing overlap for both molecules
        chem_subst_Mol1 = Chem.ReplaceCore(mol1, Chem.MolFromSmarts(overlap.smartsString), labelByIndex=True)
        chem_subst_Mol2 = Chem.ReplaceCore(mol2, Chem.MolFromSmarts(overlap.smartsString), labelByIndex=True)

        if (chem_subst_Mol1 != None and chem_subst_Mol2 != None):
            chem_substitutions_mol1 = Chem.MolToSmiles(chem_subst_Mol1)
            chem_substitutions_mol2 = Chem.MolToSmiles(chem_subst_Mol2)

            # Summarize substituents
            subst_summary = list(map(
                summarize_substitutians_across_molecule,
                [chem_substitutions_mol1, chem_substitutions_mol2]
            ))
        else:
            subst_summary = [{"": None}]

    elif overlap.numAtoms < min_overlapping_atoms:
        overlap_dict = {"smartsString": overlap.smartsString,
                        "numAtoms": overlap.numAtoms,
                        "numBonds": overlap.numBonds,
                        "canceled": overlap.canceled}

        subst_summary = [{"": None}]

    elif overlap.canceled == True:
        overlap_dict = {"smartsString": overlap.smartsString,
                        "numAtoms": overlap.numAtoms,
                        "numBonds": overlap.numBonds,
                        "canceled": overlap.canceled}

        subst_summary = [{"": None}]

    # calculate fingerprints for both SMILES
    fingerprints_classic = [Chem.RDKFingerprint(x, maxPath=fp_maxPath, fpSize=fp_fpSize) for x in [mol1, mol2]]
    fingerprints_layered = [Chem.LayeredFingerprint(x, maxPath=fp_maxPath, fpSize=fp_fpSize) for x in [mol1, mol2]]

    output = {"input": [id1, id2, smiles1, smiles2, min_overlapping_atoms],
              "overlap": overlap_dict,
              "matching_atoms": [overlap.numAtoms / mol1.GetNumAtoms() * 100,
                                 overlap.numAtoms / mol2.GetNumAtoms() * 100],
              "chemical_substituents": subst_summary,
              "similarity_metrics": {
                  "TanimotoSimilarity": [
                      DataStructs.FingerprintSimilarity(fingerprints_classic[0], fingerprints_classic[1],
                                                        metric=DataStructs.TanimotoSimilarity),
                      DataStructs.FingerprintSimilarity(fingerprints_layered[0], fingerprints_layered[1],
                                                        metric=DataStructs.TanimotoSimilarity)],
                  "DiceSimilarity": [
                      DataStructs.FingerprintSimilarity(fingerprints_classic[0], fingerprints_classic[1],
                                                        metric=DataStructs.DiceSimilarity),
                      DataStructs.FingerprintSimilarity(fingerprints_layered[0], fingerprints_layered[1],
                                                        metric=DataStructs.DiceSimilarity)],
                  "CosineSimilarity": [
                      DataStructs.FingerprintSimilarity(fingerprints_classic[0], fingerprints_classic[1],
                                                        metric=DataStructs.CosineSimilarity),
                      DataStructs.FingerprintSimilarity(fingerprints_layered[0], fingerprints_layered[1],
                                                        metric=DataStructs.CosineSimilarity)],
                  "SokalSimilarity": [
                      DataStructs.FingerprintSimilarity(fingerprints_classic[0], fingerprints_classic[1],
                                                        metric=DataStructs.SokalSimilarity),
                      DataStructs.FingerprintSimilarity(fingerprints_layered[0], fingerprints_layered[1],
                                                        metric=DataStructs.SokalSimilarity)],
                  "RusselSimilarity": [
                      DataStructs.FingerprintSimilarity(fingerprints_classic[0], fingerprints_classic[1],
                                                        metric=DataStructs.RusselSimilarity),
                      DataStructs.FingerprintSimilarity(fingerprints_layered[0], fingerprints_layered[1],
                                                        metric=DataStructs.RusselSimilarity)],
                  "KulczynskiSimilarity": [
                      DataStructs.FingerprintSimilarity(fingerprints_classic[0], fingerprints_classic[1],
                                                        metric=DataStructs.KulczynskiSimilarity),
                      DataStructs.FingerprintSimilarity(fingerprints_layered[0], fingerprints_layered[1],
                                                        metric=DataStructs.KulczynskiSimilarity)],
                  "McConnaugheySimilarity": [
                      DataStructs.FingerprintSimilarity(fingerprints_classic[0], fingerprints_classic[1],
                                                        metric=DataStructs.McConnaugheySimilarity),
                      DataStructs.FingerprintSimilarity(fingerprints_layered[0], fingerprints_layered[1],
                                                        metric=DataStructs.McConnaugheySimilarity)]
              }
              }

    return output


##example
#print(smiles_comparisson("C[CH](CCC(=O)NCC(=O)O)[CH]1CCC2[C]1([CH](CC3C2CC[CH]4[C]3(CC[CH](C4)O)C)O)C",
#                         "C[CH](CCC(=O)O)[CH]1CC[CH]2[C]1(CC[CH]3[CH]2CC[CH]4[C]3(CC[CH](C4)O)C)C",
#                         "CCMSLIB00005435519",
#                         "CCMSLIB00005435551", MCStimeout=5))
#
#print(smiles_comparisson(
#    "O[C@H](C1)CC[C@@]2(C)[C@@]1([H])C[C@@H](O)[C@]3([H])[C@]2([H])C[C@H](O)[C@@]4(C)[C@@]3([H])CC[C@@]4([C@H](C)CCCC(CO)CO)[H]",
#    "C[C@H](CCC(N[C@H](C(O)=O)CCC(O)=O)=O)[C@H]1CC[C@@]2([H])[C@]3([H])[C@H](O)C[C@]4([H])C[C@H](O)CC[C@]4(C)[C@H]3C[C@H](O)[C@@]21C",
#    "CCMSLIB00005435519",
#    "CCMSLIB00005435551", MCStimeout=5))