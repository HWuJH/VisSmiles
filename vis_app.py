import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import py3Dmol

# 设置页面标题
st.title("🔬 SMILES 处理工具")

# 侧边栏用于输入 SMILES
st.sidebar.header("📝 输入 SMILES")
smiles = st.sidebar.text_input("请输入 SMILES:", value="CCO")  # 默认乙醇 (CCO)

# 分子转换函数
def get_mol(smiles):
    """将 SMILES 转换为 RDKit 分子对象"""
    try:
        return Chem.MolFromSmiles(smiles) if smiles else None
    except:
        return None

# --- 创建布局 ---
col1, col2 = st.columns([1, 2])  # 左侧用于按钮，右侧用于展示结果

# --- 处理 SMILES ---
mol = get_mol(smiles)

with col1:  # 左侧：操作按钮
    st.subheader("📌 操作")
    
    if st.button("规范化 SMILES"):
        if mol:
            canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
            st.session_state["canonical_smiles"] = canonical_smiles  # 存储 SMILES
        else:
            st.session_state["canonical_smiles"] = "❌ 无效的 SMILES"

    if st.button("显示 2D 结构"):
        if mol:
            img = Draw.MolToImage(mol, size=(300, 300))
            st.session_state["mol_2d"] = img  # 存储 2D 结构
        else:
            st.session_state["mol_2d"] = None

    if st.button("显示 3D 结构"):
        if mol:
            try:
                mol_3d = Chem.AddHs(mol)
                if AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG()) == 0:  # 3D 坐标生成
                    mol_block = Chem.MolToMolBlock(mol_3d)
                    viewer = py3Dmol.view(width=500, height=400)
                    viewer.addModel(mol_block, "mol")
                    viewer.setStyle({"stick": {}})
                    viewer.zoomTo()
                    st.session_state["mol_3d"] = viewer._repr_html_()  # 存储 3D 结构
                else:
                    st.session_state["mol_3d"] = "⚠️ 3D 坐标生成失败"
            except Exception as e:
                st.session_state["mol_3d"] = f"⚠️ 3D 可视化失败: {e}"

# --- 右侧展示区 ---
with col2:
    st.subheader("📌 结果展示")
    
    # 规范化 SMILES 显示
    st.markdown("**✅ 规范化 SMILES:**")
    if "canonical_smiles" in st.session_state:
        st.code(st.session_state["canonical_smiles"], language="markdown")
    
    # 2D 结构展示
    st.markdown("**🧪 2D 分子结构:**")
    if "mol_2d" in st.session_state and st.session_state["mol_2d"]:
        st.image(st.session_state["mol_2d"], caption="2D 结构", use_column_width=False)
    
    # 3D 结构展示
    st.markdown("**🧩 3D 分子结构:**")
    if "mol_3d" in st.session_state:
        if isinstance(st.session_state["mol_3d"], str):  # 处理错误消息
            st.error(st.session_state["mol_3d"])
        else:
            st.components.v1.html(st.session_state["mol_3d"], height=400)
