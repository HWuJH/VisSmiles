import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import py3Dmol

# 设置页面标题
st.title("SMILES 处理工具")

# 用户输入 SMILES
smiles = st.text_input("请输入 SMILES 字符串:", value="CCO")  # 默认是乙醇 (CCO)

# 分子转换函数
def get_mol(smiles):
    """将 SMILES 转换为 RDKit 分子对象"""
    return Chem.MolFromSmiles(smiles) if smiles else None

# 规范化 SMILES
if st.button("规范化 SMILES"):
    mol = get_mol(smiles)
    if mol:
        canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
        st.success(f"规范化后的 SMILES: `{canonical_smiles}`")
    else:
        st.error("无效的 SMILES 字符串！请检查输入。")

# 2D 可视化
if st.button("显示 2D 结构"):
    mol = get_mol(smiles)
    if mol:
        img = Draw.MolToImage(mol, size=(300, 300))  # 生成 2D 图像
        st.image(img, caption="2D 分子结构", use_column_width=False)
    else:
        st.error("无效的 SMILES 字符串！请检查输入。")

# 3D 可视化
if st.button("显示 3D 结构"):
    mol = get_mol(smiles)
    if mol:
        try:
            mol_3d = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG())
            mol_block = Chem.MolToMolBlock(mol_3d)

            # 使用 py3Dmol 渲染 3D 结构
            viewer = py3Dmol.view(width=500, height=400)
            viewer.addModel(mol_block, "mol")
            viewer.setStyle({"stick": {}})
            viewer.zoomTo()
            st.components.v1.html(viewer._repr_html_(), height=400)
        except Exception as e:
            st.error(f"3D 可视化失败: {e}")
    else:
        st.error("无效的 SMILES 字符串！请检查输入。")
