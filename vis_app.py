import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import nglview as nv
from rdkit.Chem import rdmolfiles
import os

# 设置页面标题
st.title("🔬 SMILES 处理工具")

# 侧边栏输入
st.sidebar.header("📝 输入 SMILES")
smiles = st.sidebar.text_input("请输入 SMILES:", value="CCO")  # 默认乙醇 (CCO)

# 获取分子对象的函数
def get_mol(smiles):
    """将 SMILES 转换为 RDKit 分子对象"""
    try:
        return Chem.MolFromSmiles(smiles) if smiles else None
    except:
        return None

# 处理 SMILES
mol = get_mol(smiles)

# 侧边栏操作区
st.sidebar.header("📌 操作")

if st.sidebar.button("规范化 SMILES"):
    if mol:
        st.session_state["canonical_smiles"] = Chem.MolToSmiles(mol, canonical=True)
    else:
        st.session_state["canonical_smiles"] = "❌ 无效的 SMILES"

if st.sidebar.button("显示 2D 结构"):
    if mol:
        img = Draw.MolToImage(mol, size=(300, 300))
        st.session_state["mol_2d"] = img
    else:
        st.session_state["mol_2d"] = None

# 3D 渲染样式选择
render_style = st.sidebar.selectbox("选择 3D 渲染样式", ["ball+stick", "stick", "surface", "cartoon"])

if st.sidebar.button("显示 3D 结构"):
    if mol:
        try:
            mol_3d = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG())  # 生成 3D 坐标
            pdb_block = rdmolfiles.MolToPDBBlock(mol_3d)
            
            # 将 PDB 数据保存到临时文件
            with open("temp.pdb", "w") as pdb_file:
                pdb_file.write(pdb_block)
            
            # 使用 nglview 渲染 3D 分子结构
            view = nv.show_file("temp.pdb")
            
            # 设置渲染样式
            if render_style == "ball+stick":
                view.add_ball_and_stick()
            elif render_style == "stick":
                view.add_stick()
            elif render_style == "surface":
                view.add_surface()
            elif render_style == "cartoon":
                view.add_cartoon()
            
            # 将 nglview 渲染的 HTML 内容保存到 session_state
            st.session_state["mol_3d_html"] = view._repr_html_()
            
            # 清理临时文件
            os.remove("temp.pdb")
            
        except Exception as e:
            st.session_state["mol_3d_html"] = f"⚠️ 3D 可视化失败: {e}"
    else:
        st.session_state["mol_3d_html"] = None

# **调整分区布局**
col1, col2, col3 = st.columns([1.2, 1, 1.5])  # 让 3D 结构区域更大

# **规范化 SMILES 显示**
with col1:
    st.subheader("✅ 规范化 SMILES")
    if "canonical_smiles" in st.session_state:
        st.code(st.session_state["canonical_smiles"], language="markdown")

# **2D 结构显示**
with col2:
    st.subheader("🧪 2D 结构")
    if "mol_2d" in st.session_state and st.session_state["mol_2d"]:
        st.image(st.session_state["mol_2d"], caption="2D 结构", use_container_width=True)

# **3D 结构显示**
with col3:
    st.subheader("🧩 3D 结构")
    if "mol_3d_html" in st.session_state:
        if isinstance(st.session_state["mol_3d_html"], str) and "⚠️" in st.session_state["mol_3d_html"]:  # 处理错误信息
            st.error(st.session_state["mol_3d_html"])
        else:
            # 显示 nglview 渲染的 3D 分子结构
            st.components.v1.html(st.session_state["mol_3d_html"], width=500, height=400, scrolling=True)
