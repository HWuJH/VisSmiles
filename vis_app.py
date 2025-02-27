import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import py3Dmol

# 设置页面标题
st.title("🔬 SMILES 处理工具")

# 在侧边栏输入 SMILES
st.sidebar.header("📝 输入 SMILES")
smiles = st.sidebar.text_input("请输入 SMILES:", value="CCO")  # 默认乙醇 (CCO)

# 分子转换函数
def get_mol(smiles):
    """将 SMILES 转换为 RDKit 分子对象"""
    try:
        return Chem.MolFromSmiles(smiles) if smiles else None
    except:
        return None

# 规范化 SMILES
st.sidebar.subheader("📌 规范化")
if st.sidebar.button("规范化 SMILES"):
    mol = get_mol(smiles)
    if mol:
        canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
        st.success(f"✅ 规范化后的 SMILES:\n`{canonical_smiles}`")
    else:
        st.error("❌ 无效的 SMILES，请检查输入！")

# 2D 可视化
st.sidebar.subheader("🖼️ 2D 结构")
if st.sidebar.button("显示 2D 结构"):
    mol = get_mol(smiles)
    if mol:
        img = Draw.MolToImage(mol, size=(300, 300))  # 生成 2D 图像
        st.image(img, caption="🧪 2D 分子结构", use_container_width=True)
    else:
        st.error("❌ 无效的 SMILES，请检查输入！")

# 3D 可视化
st.sidebar.subheader("🧩 3D 结构")
if st.sidebar.button("显示 3D 结构"):
    mol = get_mol(smiles)
    if mol:
        try:
            mol_3d = Chem.AddHs(mol)
            if AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG()) == 0:  # 生成 3D 坐标
                mol_block = Chem.MolToMolBlock(mol_3d)

                # 3D 可视化
                viewer = py3Dmol.view(width=500, height=400)
                viewer.addModel(mol_block, "mol")
                viewer.setStyle({"stick": {}})
                viewer.zoomTo()
                st.components.v1.html(viewer._repr_html_(), height=400)
            else:
                st.error("⚠️ 3D 坐标生成失败，可能是 SMILES 过于复杂！")
        except Exception as e

