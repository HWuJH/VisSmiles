import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import py3Dmol

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

# --- 侧边栏操作区 ---
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

if st.sidebar.button("显示 3D 结构"):
    if mol:
        try:
            mol_3d = Chem.AddHs(mol)
            if AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG()) == 0:
                mol_block = Chem.MolToMolBlock(mol_3d)

                # 3D视图
                viewer = py3Dmol.view(width=500, height=400)
                viewer.addModel(mol_block, "mol")
                viewer.setStyle({"stick": {}})
                viewer.zoomTo()

                # 存储HTML代码
                st.session_state["mol_3d"] = viewer._repr_html_()
                st.session_state["mol_3d_debug"] = mol_block   # 调试信息
            else:
                st.session_state["mol_3d"] = "⚠️ 3D 坐标生成失败"
        except Exception as e:
            st.session_state["mol_3d"] = f"⚠️ 3D 可视化失败: {e}"

# --- 右侧分区调整 ---
st.markdown("<br>", unsafe_allow_html=True)    # 增加间距
col1, col2, col3 = st.columns([1, 1.2, 1.5])   # 调整列宽比例
# col1, col2, col3 = st.columns(3)

# **规范化 SMILES 显示**
with col1:
    st.subheader("✅ 规范化 SMILES")
    st.markdown("<br>", unsafe_allow_html=True)    # 增加上下间距
    if "canonical_smiles" in st.session_state:
        st.code(st.session_state["canonical_smiles"], language="markdown")

# **2D 结构显示**
with col2:
    st.subheader("🧪 2D 结构")
    st.markdown("<br>", unsafe_allow_html=True)    # 增加上下间距
    if "mol_2d" in st.session_state and st.session_state["mol_2d"]:
        st.image(st.session_state["mol_2d"], caption="2D 结构", use_container_width=True)

# **3D 结构显示**
with col3:
    st.subheader("🧩 3D 结构")
    st.markdown("<br>", unsafe_allow_html=True)    # 增加上下间距
    if "mol_3d" in st.session_state:
        if isinstance(st.session_state["mol_3d"], str):  # 处理错误信息
            st.error(st.session_state["mol_3d"])
        else:
            st.components.v1.html(st.session_state["mol_3d"], height=400, crolling=False)

    # 调试信息
    if "mol_3d_debug" in st.session_state:
            with st.expander("🛠️ 3D 结构调试信息"):
                st.text(st.session_state["mol_3d_debug"])
