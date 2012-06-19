package MiscTools.GENCODEPPCompariosonTool;


/**
 * Translation Match stores two translation objects, and is a simple abstraction of a pair of translation with the same ID that come from different sources.
 * @author David "Corvette" Thomas
 *
 */
public class TranslationMatch {

	
	private Translation gencodeTrans = null;
	private Translation ppTrans = null;
	/**
	 * @return the gencodeTrans
	 */
	public Translation getGencodeTrans() {
		return gencodeTrans;
	}
	/**
	 * @return the ppTrans
	 */
	public Translation getPpTrans() {
		return ppTrans;
	}
	/**
	 * @param gencodeTrans the gencodeTrans to set
	 */
	public void setGencodeTrans(Translation gencodeTrans) {
		this.gencodeTrans = gencodeTrans;
	}
	/**
	 * @param ppTrans the ppTrans to set
	 */
	public void setPpTrans(Translation ppTrans) {
		this.ppTrans = ppTrans;
	}
	
	
}
