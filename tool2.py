from warnings import catch_warnings
import streamlit as st
from streamlit_option_menu import option_menu
import copy
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from IPython.display import display
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG


st. set_page_config(layout="wide")
st.markdown("""
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bulma@0.9.3/css/bulma.min.css">
<nav class="navbar nav-style" role="navigation" aria-label="main navigation">
	<div class="set full navbar-brand" id="navbar">
		<a href="/" class="navbar-item false ">
			<!-- <label class="rs-phone" for="title">Powered by</label> -->
		</a>
		<img class="img2" src="https://www.mitbio.edu.in/wp-content/themes/mitsbrs/assets/images/logo.png" alt="https://www.mitbio.edu.in/wp-content/themes/mitsbrs/assets/images/logo.png"/>
	</div>
    <button class="menu-toggle">
		<span class="bar" />
		<span class="bar" />
		<span class="bar" />
	</button>
	<div class="navbar-end">
		<div class="navbar-item">
			<div class="rs-align buttons Parentbutton">
				<a href="http://localhost:3000/" onclick="func(0)" class="button bcolor false">
					<strong>Home</strong>
				</a>
				<a href="http://localhost:3000/About" class="button bcolor false">
					<strong>About</strong>
				</a>
				<a href="#sb" class="button bcolor false">
					<strong>Search</strong>
				</a>
				<a href="http://localhost:3000/tools" class="button bcolor false">
					<strong>Tools</strong>
				</a>
				<a href="http://localhost:3000/Help" class="button bcolor false">
					<strong>Help</strong>
				</a>
			</div>
		</div>
	</div>
    </nav>
""", unsafe_allow_html=True)
st.markdown("""
<style>
	.nav-style {
		background-color: rgba(0, 0, 0, 0.1);
		height: 15vh;
		padding-left: 5%;
		padding-right: 5%;
	}
	.set {
		position: relative;
		display: flex;
		flex-direction: column;
		align-items: flex-start;
	}
	.img2 {
		position: relative;
		top: -10%;
		left: 0%;
		width: 80%;
		height: 80%;
	}
	.rs-align {
		display: flex;
	}

	.menu-toggle {
		top: 0.7rem;
		right: 2%;
		position: absolute;
		display: none;
		flex-direction: column;
		justify-content: space-between;
		height: 30px;
		width: 40px;
		background: none;
		padding: 0%;
		margin: 0%;
		border: none;
	}
	.bar {
		width: 100%;
		border-radius: 10%;
		border: 2px solid rgb(36, 36, 36);
	}
	.false {
		text-decoration: none;
	}
	.bcolor {
		background: none;
		border: none;
        color:black;
	}
	.bcolor:hover {
		color: white;
	}
	.rs-phone {
		color: black;
		font-weight: 700;
		font-size: 8px;
	}

	@media screen and (max-width: 1050px) {
		.set {
			display: flex;
			flex-direction: column;
			align-items: flex-start;
			height: 100%;
		}
		.rs-phone {
			font-size: 8px;
		}
		.menu-toggle {
			display: flex;
		}
		.rs-align {
			display: none;
		}
		.img2 {
			height: 100px;
			width: 200px;
		}

	}
	@media screen and (max-width:1100px){
		.nav-style{
			padding-left: 0%;
			padding-right: 0%;
		}
	}
</style>""", unsafe_allow_html=True)

st.markdown("""<style>
img{

}
</style>""", unsafe_allow_html=True)

st.markdown("""
        <style>
               .css-18e3th9 {
                    padding-top: 0rem;
                    padding-bottom: 10rem;
                    padding-left: 0rem;
                    padding-right: 0rem;
                }
               .css-1d391kg {
                    padding-top: 3.5rem;
                    padding-right: 1rem;
                    padding-bottom: 3.5rem;
                    padding-left: 1rem;
                }
        </style>
        """, unsafe_allow_html=True)

hide_st_style = """
            <style>
            # MainMenu {visibility: hidden;}
            footer {visibility: hidden;}
            header {visibility: hidden;}
            </style>
            """
st.markdown(hide_st_style, unsafe_allow_html=True)

smiles = st.text_input('Enter SMILES below:',
                       'CC(C)(C)C1=CC=C(C=C1)C(=O)CC(=O)C2=CC=C(C=C2)OC')
st.markdown(
    """
<style>
.row-widget{
position:reletive;
font-weight:600;
left:2%;
}
</style>
""", unsafe_allow_html=True)
# processing starts
try:
    smi = smiles
    mol_t = Chem.MolFromSmiles(smi)
    im1 = Draw.MolToImage(mol_t)
    st.image(im1)
except:
    st.markdown(
        """
<h1 class="valid">Please enter a valid SMILES String</h1>
""", unsafe_allow_html=True)
st.markdown(
    """
<style>
.valid{
position:reletive;
left:2%;
}
</style>

""", unsafe_allow_html=True)
for i in mol_t.GetAtoms():
    i.SetIntProp("atom_idx", i.GetIdx())
for i in mol_t.GetBonds():
    i.SetIntProp("bond_idx", i.GetIdx())
mol_s = MurckoScaffold.GetScaffoldForMol(mol_t)

im2 = Draw.MolToImage(mol_s)
st.image(im2)

ri = mol_s.GetRingInfo()
   # You can interrogate the RingInfo object to tell you the atoms that make up each ring:
print(ri.AtomRings())

ri = mol_t.GetRingInfo()
bondrings = ri.BondRings()

bondring_list = list(bondrings[0]+bondrings[1])


all_bonds_idx = [bond.GetIdx() for bond in mol_t.GetBonds()]


none_ring_bonds_list = []
for i in all_bonds_idx:
        if i not in bondring_list:
            none_ring_bonds_list.append(i)

cut_bonds = []
for bond_idx in none_ring_bonds_list:
        bgn_atom_idx = mol_t.GetBondWithIdx(bond_idx).GetBeginAtomIdx()
        ebd_atom_idx = mol_t.GetBondWithIdx(bond_idx).GetEndAtomIdx()
        if mol_t.GetBondWithIdx(bond_idx).GetBondTypeAsDouble() == 1.0:
            if mol_t.GetAtomWithIdx(bgn_atom_idx).IsInRing()+mol_t.GetAtomWithIdx(ebd_atom_idx).IsInRing() == 1:
                t_bond = mol_t.GetBondWithIdx(bond_idx)
                t_bond_idx = t_bond.GetIntProp("bond_idx")
                cut_bonds.append(t_bond_idx)
res = Chem.FragmentOnBonds(mol_t, cut_bonds)
im4 = Draw.MolToImage(res)

st.image(im4)

frgs = Chem.GetMolFrags(res, asMols=True)
im5 = Draw.MolsToImage(frgs)
st.image(im5)


# footer
st.markdown(
   
   '''
<div style="background:black" >
        <div class="content2" style="display: flex; justify-content: center;
  align-items: center;">
		<p style="margin-top:2.5%; color:white"
		>
			Copyright &copy MIT ADT, All rights reserved.
		</p>
	</div>
</div>
   '''
   , unsafe_allow_html=True)

st.markdown("""<style>
    .css-qri22k{
        padding-top: 0rem;
        padding-bottom:0rem;
        padding-left: 0rem;
        padding-right: 0rem;}
        </style> """, unsafe_allow_html=True) 
